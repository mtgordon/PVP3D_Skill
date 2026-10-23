"""Batch 286 (2026-10-15; the user, on L283_junction_fine2z18_MAX's curves folding through themselves: "could we remesh the
curve so that the thickness stays the same, but the other two dimensions are meshed more closely together?"; picture
curve_refine_view/refine_around_right.png approved: "build it as shown and run MAX").
NOT IN SOURCE (mesh refinement only): every curve hex (canal_wrap_hex_avw / _pvw / _pvw_fine) with its centre at rest
z < ZCUT is split in 2 AROUND the curve's U: its edges 0-1, 3-2, 4-5, 7-6 (~2 mm; 1-2 = through the wall, 4 layers of
~0.6 mm, untouched; 0-4 = along the canal, already halved by build_batch280) get a mid node. The first half keeps the
element id and domain, the second half goes into the new domain canal_wrap_hex_around (Curve_soft, same mass damping).
New nodes: on an edge also used by an element that is NOT split (the zone border at ZCUT, or a wall element sharing the
curve's nodes) -> held at the edge midpoint (linear constraints, curve_around_ties); otherwise on an edge whose two ends are
on / tied onto the PVW, PeB or AVW -> tied onto that part's face (as build_batch280; junction_ties); otherwise free.
Surface facets on split hex faces are split in 2 (every surface: contacts, the curve pressure, ...).
usage: py -3.10 build_batch286.py SRC_RUN NEW_RUN [--zcut -15]"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch278 import bilinear, lc, PVW, PEB, EDGES  # noqa: E402
from apex_close import HF  # noqa: E402
from build_batch162 import write  # noqa: E402
from curve_swell import TETS  # noqa: E402

DOMS = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
AROUND = ((0, 1), (3, 2), (4, 5), (7, 6))


def main(src, new, zcut=-15.0):
    t, X, E = load(src)
    doms = DOMS + tuple(sorted({n for n, _ in E.values() if n.startswith('canal_wrap_hex') and n.endswith('_len')
                                 and n != 'canal_wrap_hex_avw_cervix_len'}))   # 10-23: also build_batch412's along-canal halves
    zone = [k for k, (n, c) in E.items() if n in doms and np.mean([X[v][2] for v in c]) < zcut]
    S = {frozenset((E[k][1][a], E[k][1][b])) for k in zone for a, b in AROUND}
    # every element's edges (hex8: the 12 edges; others: consecutive nodes)
    zs = set(zone)
    users = defaultdict(set)
    for k, (n, c) in E.items():
        ed = [(c[a], c[b]) for a, b in EDGES] if len(c) == 8 else [(c[i], c[(i + 1) % len(c)]) for i in range(len(c))]
        for a, b in ed:
            e = frozenset((a, b))
            if e in S:
                users[e].add(k)
    # consistency: a split edge must be an AROUND edge of every zone hex that uses it
    for k in zone:
        c = E[k][1]
        ar = {frozenset((c[a], c[b])) for a, b in AROUND}
        bad = [e for e in (frozenset((c[a], c[b])) for a, b in EDGES) if e in S and e not in ar]
        assert not bad, f'hex {k}: a split edge is not one of its around edges (mixed orientation)'
    hanging = {e for e in S if users[e] - zs}
    # parts the nodes are on / tied onto (build_batch280's rule)
    dn = lambda d: {v for n, c in E.values() if n == d for v in c}   # noqa: E731
    SETS = {'PVW': dn(PVW), 'PeB': dn(PEB), 'AVW': dn('_PickedSet347')}
    held = {}
    for name in ('junction_ties', 'canal_wrap_ends'):
        blk = re.search(r'<constraint name="%s".*?</constraint>' % name, t, re.S)
        for lcs in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', blk.group(0), re.S):
            nn = [int(a) for a in re.findall(r'<node id="(\d+)" bc="x">', lcs)]
            if nn:
                held[nn[0]] = next((p for p, S_ in SETS.items() if set(nn[1:]) <= S_), None)
    on = lambda v, part: v in SETS[part] or held.get(v) == part   # noqa: E731
    pfaces = {}
    for part, dom in (('PVW', PVW), ('PeB', PEB), ('AVW', '_PickedSet347')):
        els = [c for n, c in E.values() if n == dom]
        cnt = defaultdict(int)
        for c in els:
            for fc in HF:
                cnt[frozenset(c[i] for i in fc)] += 1
        pfaces[part] = [[c[i] for i in fc] for c in els for fc in HF if cnt[frozenset(c[i] for i in fc)] == 1
                        or (part == 'PeB' and set(c[i] for i in fc) <= (SETS['PVW'] & SETS['PeB']))]
    nid = max(X) + 1
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    mid, pos = {}, {}
    for e in sorted(S, key=lambda e: tuple(sorted(e))):
        a, b = sorted(e)
        mid[e] = nid; pos[nid] = (X[a] + X[b]) / 2; nid += 1
    ties_mid, ties_face, nt, worst = '', '', defaultdict(int), 0.0
    for e, v in mid.items():
        a, b = sorted(e)
        if e in hanging:
            ties_mid += lc(v, [(a, 0.5), (b, 0.5)]); nt['midpoint (hanging)'] += 1
            continue
        for part in ('PVW', 'PeB', 'AVW'):
            if on(a, part) and on(b, part):
                F = pfaces[part]
                C = np.array([np.mean([X[q] for q in f], 0) for f in F])
                best = None
                for i in np.argsort(np.linalg.norm(C - pos[v], axis=1))[:30]:
                    Q = np.array([X[q] for q in F[i]])
                    w, d, m = bilinear(Q, pos[v])
                    if m <= 1.0 + 1e-3 and (best is None or d < best[1]):
                        best = (list(zip(F[i], w)), d, w @ Q)
                if best is None:
                    continue
                worst = max(worst, best[1]); pos[v] = best[2]
                ties_face += lc(v, best[0]); nt[f'onto {part} face'] += 1
                break
    # sub-hexes
    M = lambda c, a, b: mid[frozenset((c[a], c[b]))]   # noqa: E731
    sub, second = {}, []
    for k in zone:
        c = E[k][1]
        m01, m32, m45, m76 = M(c, 0, 1), M(c, 3, 2), M(c, 4, 5), M(c, 7, 6)
        A = [c[0], m01, m32, c[3], c[4], m45, m76, c[7]]
        B = [m01, c[1], c[2], m32, m45, c[5], c[6], m76]
        sub[k] = A; second.append((eid, B)); eid += 1
    P = lambda v: pos.get(v, X.get(v))   # noqa: E731
    mn = min(np.dot(P(c[q1]) - P(c[q0]), np.cross(P(c[q2]) - P(c[q0]), P(c[q3]) - P(c[q0]))) / 6
             for c in list(sub.values()) + [b for _, b in second] for q0, q1, q2, q3 in TETS)
    assert mn > 0, mn
    # ---- write: nodes (after the last Nodes block), elements, domain, damping
    b = t.rindex('</Nodes>\n') + len('</Nodes>\n')
    t = t[:b] + '\t\t<Nodes name="curve_around_nodes">\n' + ''.join(
        f'\t\t\t<node id="{v}">' + ','.join(f'{x:.7g}' for x in pos[v]) + '</node>\n' for v in sorted(pos)) + '\t\t</Nodes>\n' + t[b:]
    for k in zone:
        old = f'<elem id="{k}">' + ','.join(map(str, E[k][1])) + '</elem>'
        assert old in t, k
        t = t.replace(old, f'<elem id="{k}">' + ','.join(map(str, sub[k])) + '</elem>', 1)
    a = t.index('<Elements type="hex8" name="canal_wrap_hex_pvw_fine">'); b = t.index('</Elements>', a) + len('</Elements>\n')
    t = t[:b] + '\t\t<Elements type="hex8" name="canal_wrap_hex_around">\n' + ''.join(
        f'\t\t\t<elem id="{i}">' + ','.join(map(str, c)) + '</elem>\n' for i, c in second) + '\t\t</Elements>\n' + t[b:]
    t = t.replace('\t\t<SolidDomain name="canal_wrap_hex_pvw_fine" mat="Curve_soft" />\n',
                  '\t\t<SolidDomain name="canal_wrap_hex_pvw_fine" mat="Curve_soft" />\n'
                  '\t\t<SolidDomain name="canal_wrap_hex_around" mat="Curve_soft" />\n', 1)
    m = re.search(r'\t\t<body_load name="wall_damping_canal_wrap_hex_pvw_fine" .*?</body_load>\n', t, re.S)
    t = t[:m.end()] + m.group(0).replace('canal_wrap_hex_pvw_fine', 'canal_wrap_hex_around') + t[m.end():]
    # surfaces: a quad facet with two split edges on opposite sides -> 2 quads, same orientation
    nsplit = 0

    def surf(mm):
        nonlocal nsplit
        out, n = [], 0
        for fm in re.finditer(r'<(quad4|tri3) id="\d+">([\d,]+)</\1>', mm.group(2)):
            f = [int(v) for v in fm.group(2).split(',')]
            sp = [i for i in range(len(f)) if frozenset((f[i], f[(i + 1) % len(f)])) in mid] if fm.group(1) == 'quad4' else []
            if sp in ([0, 2], [1, 3]):
                r = sp[0]; p = f[r:] + f[:r]                    # p0-p1 and p2-p3 split
                m1, m2 = mid[frozenset((p[0], p[1]))], mid[frozenset((p[2], p[3]))]
                for g in ([p[0], m1, m2, p[3]], [m1, p[1], p[2], m2]):
                    n += 1; out.append(f'\t\t\t<quad4 id="{n}">' + ','.join(map(str, g)) + '</quad4>\n')
                nsplit += 1
            else:
                # a face of an unsplit element on the zone border: its split edge is held at the midpoint
                assert all(frozenset((f[i], f[(i + 1) % len(f)])) in hanging for i in sp), (mm.group(1), f, sp)
                n += 1; out.append(f'\t\t\t<{fm.group(1)} id="{n}">{fm.group(2)}</{fm.group(1)}>\n')
        return mm.group(1) + ''.join(out) + '\t\t</Surface>'
    t = re.sub(r'(<Surface name="[^"]+">\n)(.*?)\t\t</Surface>', surf, t, flags=re.S)
    blk = re.search(r'\t\t<constraint name="junction_ties".*?</constraint>\n', t, re.S).group(0)
    head = re.match(r'(\t\t<constraint name="junction_ties".*?\n)(.*?)\t\t\t<linear_constraint>', blk, re.S)
    around = (head.group(1).replace('junction_ties', 'curve_around_ties') + head.group(2) + ties_mid + '\t\t</constraint>\n'
              if ties_mid else '')
    t = t.replace(blk, blk.replace('\t\t</constraint>\n', ties_face + '\t\t</constraint>\n') + around, 1)
    # other constraints whose hosts are nodes of a split face (their bilinear weights now span a split face)
    hostnote = sum(1 for lcs in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', t, re.S)
                   if len(set(int(a) for a in re.findall(r'<node id="(\d+)" bc="x">', lcs)[1:]) & {v for e in S for v in e}) >= 2)
    msg = (f'NOT IN SOURCE (mesh refinement): every curve hex ({", ".join(DOMS)}) with its centre at rest z < {zcut} '
           f'({len(zone)} hexes) split in 2 around the curve\'s U (through-wall layers unchanged) -> {2 * len(zone)} hexes '
           f'({len(second)} new, domain canal_wrap_hex_around, Curve_soft, mass damping as the curves); {len(mid)} new nodes: '
           f'{dict(nt)} (face ties moved <= {worst:.3f} mm), the rest free; {nsplit} surface facets split; '
           f'{hostnote} existing tie rows have >= 2 hosts on split edges (kept as they were); smallest sub-tet {mn:.4f} mm3')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch286.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    zc = float(a[a.index('--zcut') + 1]) if '--zcut' in a else -15.0
    main(a[0], a[1], zc)
