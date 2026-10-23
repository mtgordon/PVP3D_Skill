"""Batch 403 (2026-10-22; the user: "it looks like the straps were extruded back from the PM arch. If they were extruded forward
... could you then attach them to the distal end of the curves?" -> yes). Two changes to the fan strap at the curve tip, both
NOT IN SOURCE (the strap and the curves are ours):
 (1) FLIP: the strap (PM_avw_bottom_left/right_sheet_solid, one penta6 layer TH ~0.49 mm) was extruded from its base sheet back
     INTO the curve (its first wedges overlap the curve's end hexes: the knife-edge corner). Its top nodes are moved to the other
     side of the base (top' = base - (top - base)) and each wedge is renumbered (c3, c4, c5, c0, c1, c2) so it stays positive.
     The base nodes (the curve-rim tie row, the junction / patch ties, the PM-end anchors) do not move; element ids, fibre
     directions (mat_axis) and materials are kept; every surface facet made only of strap nodes is re-wound (the flip reverses
     their outward normal).
 (2) CAP: every free distal end face of the curves (boundary quads within 1.5 mm of the end along the end tangent, normal along
     it; faces with a node shared with a domain other than the curves / strap, i.e. the PVW / PeB junction, left out) is
     extruded distally by TH into one hex8 layer (domain strap_cap, material strap_cap_mat = the strap's matrix: Mooney-Rivlin
     c1 C1, c2 0, k K, no fibre) sharing the end-face nodes (tied over the whole face instead of one rim row). Where an end-face
     node is a strap base node, the cap uses that node's strap top partner, so cap and strap join along the rim. Every surface
     facet equal to a capped end face is replaced by the cap's outer face (same winding); the cap's free side faces are added
     to curve_all (the curve self-contact).
usage: py -3.10 build_batch403.py SRC_RUN NEW_RUN [--th 0.49] [--dry]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402

CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
STRAP = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]
W3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def tv(P, c, split):
    return sum(np.dot(np.cross(P[c[b]] - P[c[a]], P[c[d]] - P[c[a]]), P[c[e]] - P[c[a]]) / 6 for a, b, d, e in split)


def main(src, new, th, dry, src_r=0.3):
    t, X, E = load(src)
    X = {k: v.copy() for k, v in X.items()}
    dn = defaultdict(set)
    for nm, c in E.values():
        for v in c:
            dn[v].add(nm)
    # ---- (1) flip the strap
    wed = {k: c for k, (nm, c) in E.items() if nm in STRAP}
    partner = {}
    for c in wed.values():
        for i in range(3):
            assert partner.get(c[i], c[i + 3]) == c[i + 3]
            partner[c[i]] = c[i + 3]
    base = set(partner); top = set(partner.values())
    assert not base & top
    assert all(dn[v] <= set(STRAP) for v in top), 'a strap top node is shared with another domain'
    for b_, t_ in partner.items():
        X[t_] = 2 * X[b_] - X[t_]
    newwed = {k: [c[3], c[4], c[5], c[0], c[1], c[2]] for k, c in wed.items()}
    bad = [k for k, c in newwed.items() if tv(X, c, W3) <= 0]
    # ---- (2) the cap
    hexes = {k: c for k, (nm, c) in E.items() if nm in CURVES and len(c) == 8}
    cnt = Counter(frozenset(c[i] for i in f) for c in hexes.values() for f in FACES)
    caps, skipped = [], 0
    for sg in (-1, 1):
        cn = np.array([X[v] for c in hexes.values() for v in c if np.sign(X[v][0]) == sg])
        lo = cn[:, 2].min(); e = cn[cn[:, 2] < lo + 3].mean(0); u = cn[(cn[:, 2] > lo + 8) & (cn[:, 2] < lo + 12)].mean(0)
        tan = (e - u) / np.linalg.norm(e - u)
        smax = max(X[v] @ tan for c in hexes.values() for v in c if np.sign(X[v][0]) == sg)
        for c in hexes.values():
            for f in FACES:
                q = [c[i] for i in f]
                if cnt[frozenset(q)] != 1 or np.sign(np.mean([X[v][0] for v in q])) != sg:
                    continue
                if any(smax - X[v] @ tan > 1.5 for v in q):
                    continue
                P = np.array([X[v] for v in q]); n = np.cross(P[2] - P[0], P[3] - P[1]); n /= np.linalg.norm(n)
                if n @ tan < 0.8:
                    continue
                if any(dn[v] - set(CURVES) - set(STRAP) for v in q):
                    skipped += 1; continue
                caps.append((q, n))
    # cap top nodes: averaged end-face normal, or the strap top partner on the rim
    nsum = defaultdict(lambda: np.zeros(3))
    for q, n in caps:
        for v in q:
            nsum[v] += n
    nid = max(X) + 1; tnode = {}; newnodes = []
    for v in sorted(nsum):
        if v in partner:
            tnode[v] = partner[v]
        else:
            tnode[v] = nid; X[nid] = X[v] + th * nsum[v] / np.linalg.norm(nsum[v]); newnodes.append(nid); nid += 1
    caphex = [q + [tnode[v] for v in q] for q, n in caps]
    # (L402 check) a moved strap top or a new cap node inside another solid at rest (the flipped strap / cap next to the PeB):
    # shorten its offset from its base node until it is outside, with a 0.05 mm margin (the layer is locally thinner there)
    baseof = {t_: b_ for b_, t_ in partner.items()}
    baseof.update({tnode[v]: v for v in nsum if v not in partner})
    tets, tel = [], []
    for k, (nm, c) in E.items():
        if nm in STRAP:
            continue
        sp = T5 if len(c) == 8 else (W3 if len(c) == 6 else None)
        if sp:
            tets += [[c[i] for i in s] for s in sp]
    from scipy.spatial import cKDTree
    TP = np.array([[X[v] for v in q] for q in tets]); cen = TP.mean(1)
    tr = cKDTree(cen); RM = np.max(np.linalg.norm(TP - cen[:, None], axis=2))

    def inside(p, skip):
        for j in tr.query_ball_point(p, RM):
            if skip & set(tets[j]):
                continue
            a_, b_, c_, d_ = TP[j]
            lam = np.linalg.solve(np.column_stack([b_ - a_, c_ - a_, d_ - a_]), p - a_)
            if min(1 - lam.sum(), lam.min()) > 1e-6:
                return True
        return False
    shortened = []
    for v, b_ in baseof.items():
        if not inside(X[v], {v, b_}):
            continue
        off = X[v] - X[b_]; lo_, hi_ = 0.0, 1.0
        for _ in range(30):
            s_ = (lo_ + hi_) / 2
            lo_, hi_ = (s_, hi_) if not inside(X[b_] + s_ * off, {v, b_}) else (lo_, s_)
        s_ = max(lo_ - 0.05 / np.linalg.norm(off), 0.2)
        X[v] = X[b_] + s_ * off; shortened.append((v, round(s_, 2)))
    print('nodes shortened to stay outside other solids (node, fraction of the thickness kept):', shortened)
    bad = [k for k, c in newwed.items() if tv(X, c, W3) <= 0]
    badc = [i for i, c in enumerate(caphex) if tv(X, c, T5) <= 0]
    print(f'strap wedges {len(wed)} (flipped; non-positive after renumbering: {len(bad)}); cap: {len(caps)} end faces '
          f'({skipped} touching the PVW / PeB junction left out), {len(newnodes)} new nodes, rim nodes on strap tops '
          f'{sum(1 for v in nsum if v in partner)}; non-positive cap hexes {len(badc)}')
    if dry or bad or badc:
        return
    # ---- write: node coordinates (strap tops), new nodes, elements, domain, material
    def node_sub(m):
        k = int(m.group(1))
        return f'{m.group(0)[:m.start(2) - m.start(0)]}{X[k][0]:.7g},{X[k][1]:.7g},{X[k][2]:.7g}</node>' if k in top else m.group(0)
    shortened_note = f'; {len(shortened)} nodes next to the PeB kept outside it (offset scaled to {min([s for _, s in shortened] or [1]):.2f}-1 of TH)'
    t = re.sub(r'<node id="(\d+)">([^<]+)</node>', node_sub, t)
    i = t.rindex('\t\t</Nodes>\n') + len('\t\t</Nodes>\n')
    t = t[:i] + ('\t\t<Nodes name="strap_cap_nodes">\n' + ''.join(
        f'\t\t\t<node id="{v}">{X[v][0]:.7g},{X[v][1]:.7g},{X[v][2]:.7g}</node>\n' for v in newnodes) + '\t\t</Nodes>\n') + t[i:]

    def el_sub(m):
        k = int(m.group(1))
        return f'\t\t\t<elem id="{k}">{",".join(map(str, newwed[k]))}</elem>\n' if k in newwed else m.group(0)
    t = re.sub(r'\t\t\t<elem id="(\d+)">([0-9,]+)</elem>\n', el_sub, t)
    eid = max(E) + 1
    i = t.rindex('\t\t</Elements>\n') + len('\t\t</Elements>\n')
    t = t[:i] + ('\t\t<Elements type="hex8" name="strap_cap">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, c))}</elem>\n' for j, c in enumerate(caphex)) + '\t\t</Elements>\n') + t[i:]
    a = t.index('<material id="', t.index('name="PM_avw_bottom_left_sheet_mat"') - 40)
    blk = t[a:t.index('</material>', a)]
    c1 = float(re.search(r'<c1>([^<]+)</c1>', blk).group(1)); k_ = float(re.search(r'<k>([^<]+)</k>', blk).group(1))
    rho = re.search(r'<density>([^<]+)</density>', blk).group(1)
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    i = t.index('</Material>')
    t = t[:i] + (f'\t\t<material id="{mid}" name="strap_cap_mat" type="Mooney-Rivlin">\n\t\t\t<density>{rho}</density>\n'
                 f'\t\t\t<c1>{c1:g}</c1>\n\t\t\t<c2>0</c2>\n\t\t\t<k>{k_:g}</k>\n\t\t</material>\n') + t[i:]
    i = t.index('\t</MeshDomains>')
    t = t[:i] + '\t\t<SolidDomain name="strap_cap" mat="strap_cap_mat"/>\n' + t[i:]
    # ---- surfaces: re-wind all-strap facets, replace capped end faces, add cap side faces to curve_all
    strapn = base | top
    capmap = {frozenset(q): [tnode[v] for v in q] for q, n in caps}
    allf = Counter(frozenset(c[i] for i in f) for c in caphex for f in FACES)
    for c in newwed.values():
        for f in ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)):
            allf[frozenset(c[i] for i in f)] += 1
    capside = [[c[i] for i in f] for c in caphex for f in FACES[2:] if allf[frozenset(c[i] for i in f)] == 1]
    stats = Counter()
    curve_all_after = []

    def surf_sub(ms):
        name, body = ms.group(1), ms.group(2)
        out = []
        for g, fid, nodes in re.findall(r'\t\t\t<(quad4|tri3) id="(\d+)">([0-9,]+)</\1>\n', body):
            q = [int(v) for v in nodes.split(',')]
            if set(q) <= strapn:
                q = q[::-1]; stats[name, 'rewound'] += 1
            elif frozenset(q) in capmap:
                if name == 'curve_all':   # (L402 check) the 0.49 mm cap is thinner than curve_all_self's 0.5 mm radius:
                    stats[name, 'end faces removed (capped)'] += 1; continue   # the cap gets its own contact (cap_contact)
                q = capmap[frozenset(q)]; stats[name, 'capped'] += 1
            out.append((g, q))
        if name == 'curve_all':
            curve_all_after.extend(q for g, q in out)
        return (f'\t\t<Surface name="{name}">\n' + ''.join(f'\t\t\t<{g} id="{i}">{",".join(map(str, q))}</{g}>\n'
                                                       for i, (g, q) in enumerate(out, 1)) + '\t\t</Surface>\n')
    t = re.sub(r'\t\t<Surface name="([^"]+)">\n(.*?)\t\t</Surface>\n', surf_sub, t, flags=re.S)
    # (L402 check) curve nodes on the cap: a separate contact as L399c (build_batch399c): primary = curve_all faces with no
    # cap / strap node, secondary = the cap's outer faces (tops + free sides), search radius SRC < the cap thickness
    capn = {v for c in caphex for v in c}
    prim = [q for q in curve_all_after if not set(q) & (capn | strapn)]
    sec = list(capmap.values()) + capside
    gq = lambda f: 'tri3' if len(f) == 3 else 'quad4'   # noqa: E731
    surf = ('\t\t<Surface name="cap_contact_primary">\n' + ''.join(
        f'\t\t\t<{gq(f)} id="{k}">{",".join(map(str, f))}</{gq(f)}>\n' for k, f in enumerate(prim, 1)) +
        '\t\t</Surface>\n\t\t<Surface name="cap_contact_secondary">\n' + ''.join(
        f'\t\t\t<quad4 id="{k}">{",".join(map(str, f))}</quad4>\n' for k, f in enumerate(sec, 1)) +
        '\t\t</Surface>\n\t\t<SurfacePair name="cap_contact">\n\t\t\t<primary>cap_contact_primary</primary>\n'
        '\t\t\t<secondary>cap_contact_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    a = t.index('<contact name="curve_all_self"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = t[a:b].replace('name="curve_all_self" surface_pair="curve_all_self"', 'name="cap_contact" surface_pair="cap_contact"')
    c = re.sub(r'<search_radius>[^<]*</search_radius>', f'<search_radius>{src_r:g}</search_radius>', c)
    t = t[:b] + '\t\t' + c + t[b:]
    stats['cap_contact', f'primary {len(prim)} / secondary {len(sec)}, search_radius {src_r:g}'] = 1
    st = '; '.join(f'{s} {k} {v}' for (s, k), v in sorted(stats.items()))
    msg = (f'NOT IN SOURCE (the user 10-22): (1) the fan strap flipped to the distal side of its base sheet (top nodes mirrored '
           f'through the base, {len(wed)} wedges renumbered; base / ties / anchors unchanged); (2) strap_cap: {len(caphex)} hex8 '
           f'({th:g} mm, strap matrix Mooney-Rivlin c1 {c1:g} k {k_:g}, no fibre) over the curves\' free distal end faces '
           f'({skipped} junction faces left out), {len(newnodes)} new nodes, joined to the strap tops on the rim{shortened_note}; '
           f'surfaces: {st}')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch403.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--th') + 1]) if '--th' in a else 0.49, '--dry' in a)
