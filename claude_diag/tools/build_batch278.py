"""Batch 278 (2026-10-10, the user: "My guess would be to get rid of the flares and make the curves attach better to the PeB ...
try some things out and see what works"). The tube main's distal junction, rebuilt (NOT IN SOURCE):
  * NO FLARE: the PVW's last layer (the 78 hexes on the PVW / PeB face; build_batch136 flared it from |x| 6 at depth 5.3 to
    the PeB's full 8.6 at the face) goes straight down: each of its 108 bottom nodes, shared with the PeB until now, becomes a
    new PVW-only node straight below its top partner (lateral position of the top node, AP / depth of the face row), and is
    tied onto the PeB's top facet it lies on (exact facet ties, linear constraints, as canal_wrap_ends). The PeB keeps its
    full width and its nodes. The PVW-side surfaces (lumen contact, pressure, walls_LA) follow the new nodes.
  * CURVES ON THE PVW: the curves' last two PVW-side rings (built beside the narrow PVW, cut by the flare until now) join the
    straight PVW face as higher up: the 4 end-column nodes ARE the PVW's lateral bottom nodes (merged), the 5th and the middle
    ring's face column are tied onto the PVW's lateral face (facet / edge ties).
  * CURVES ON THE PeB (mode full): every other end-ring node of the curves that lies over the PeB's top face (inside a
    facet, within DTOL) is snapped onto it and tied to it, so the curve end sits on the PeB.
  * canal_wrap_embed (the 20 point ties inside PVW / PeB hexes) removed.
modes: full | nopeb (no curve-end ties to the PeB: isolates them)
usage: py -3.10 build_batch278.py MODE CASE [CASE ...]   (CASE: healthy, MAX, F8A, P1 -> L258_tubecand_CASE)"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load, frame, overlap  # noqa: E402
from junction_proposal import ring_faces, mirror_ids  # noqa: E402
from apex_close import HF  # noqa: E402
from build_batch162 import write  # noqa: E402

PVW, PEB = '_PickedSet64', '_PickedSet66'
RIGHT = [(27795, 27806), (27807, 27818), (27819, 27830)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
DTOL = 1.0          # mm: a curve end node this close to a PeB top facet (inside it) is snapped and tied
NAMES = {'full': 'L278_junction', 'nopeb': 'L279_junction_nopeb'}


def bilinear(P, p):
    """(weights, distance, max |xi|) of the point on quad P (4x3) nearest p (Newton on the bilinear map)"""
    xi = np.zeros(2)
    for _ in range(30):
        r, s = xi
        N = 0.25 * np.array([(1 - r) * (1 - s), (1 + r) * (1 - s), (1 + r) * (1 + s), (1 - r) * (1 + s)])
        dr = 0.25 * np.array([-(1 - s), (1 - s), (1 + s), -(1 + s)])
        ds = 0.25 * np.array([-(1 - r), -(1 + r), (1 + r), (1 - r)])
        x = N @ P; J = np.array([dr @ P, ds @ P]).T
        step = np.linalg.lstsq(J, p - x, rcond=None)[0]
        xi = xi + step
        if np.abs(step).max() < 1e-10:
            break
    r, s = xi
    N = 0.25 * np.array([(1 - r) * (1 - s), (1 + r) * (1 - s), (1 + r) * (1 + s), (1 - r) * (1 + s)])
    return N, float(np.linalg.norm(N @ P - p)), float(np.abs(xi).max())


def lc(node, hosts):
    """three linear constraints (x, y, z): node = sum w_q host_q"""
    out = ''
    for dof in 'xyz':
        out += (f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{node}" bc="{dof}">1</node>\n'
                + ''.join(f'\t\t\t\t<node id="{q}" bc="{dof}">{-w:.6f}</node>\n' for q, w in hosts if abs(w) > 1e-7)
                + '\t\t\t</linear_constraint>\n')
    return out


def build(mode, case):
    src = f'L258_tubecand_{case}'
    t, X, E = load(src)
    c0, B = frame(X, E)
    L = {k: (p - c0) @ B for k, p in X.items()}
    dn = lambda d: {v for n, c in E.values() if n == d for v in c}   # noqa: E731
    pvw, peb = dn(PVW), dn(PEB)
    S = pvw & peb
    PL = {k: c for k, (n, c) in E.items() if n == PVW and set(c) & S}
    # top partner of each shared node: its vertical hex edge in the PVW's last layer
    partner = defaultdict(set)
    for c in PL.values():
        for a, b in EDGES:
            if (c[a] in S) != (c[b] in S):
                s, u = (c[a], c[b]) if c[a] in S else (c[b], c[a])
                partner[s].add(u)
    assert all(len(v) == 1 for v in partner.values()) and len(partner) == len(S), 'last-layer topology'
    partner = {s: next(iter(v)) for s, v in partner.items()}
    newpos = {s: np.array([L[partner[s]][0], L[s][1], L[s][2]]) for s in S}
    # the curves' rings, both sides
    sides = [RIGHT, []]
    for a, b in RIGHT:
        m = mirror_ids(E, L, range(a, b + 1)); sides[1].append((min(m), max(m)))
    merged, face_ties, end_rings, log = {}, [], [], []
    for rings in sides:
        _, botA = ring_faces(E, L, rings[0])
        _, mid = ring_faces(E, L, rings[1])
        _, end = ring_faces(E, L, rings[2])
        end_rings.append(sorted(end))
        face_top = sorted(botA & pvw, key=lambda v: L[v][1])
        corner = {u: s for s, u in partner.items()}                        # top node -> its shared bottom node
        bottom = [corner[u] for u in face_top]                              # the PVW's lateral bottom column (4 rows)
        end_face = sorted(sorted(end, key=lambda v: abs(L[v][0]))[:5], key=lambda v: L[v][1])
        mid_face = sorted(sorted(mid, key=lambda v: abs(L[v][0]))[:5], key=lambda v: L[v][1])
        for s in bottom:                                                    # merge: nearest end-face node by AP
            v = min(end_face, key=lambda q: abs(L[q][1] - L[s][1]))
            assert abs(L[v][1] - L[s][1]) < 0.1, (v, s)
            merged[s] = v
        log.append(f'rings {rings}: PVW face top {face_top}, lateral bottom (shared) {bottom} -> curve nodes '
                   f'{[merged[s] for s in bottom]} (moved {np.round([np.linalg.norm(newpos[s] - L[merged[s]]) for s in bottom], 2).tolist()} mm)')
        # the PVW's straight lateral face of this side: quads (top row u1, u2 + their new bottoms)
        for v in [q for q in end_face if q not in merged.values()] + mid_face:
            best = None
            for j in range(len(face_top) - 1):
                u1, u2 = face_top[j], face_top[j + 1]
                b1, b2 = bottom[j], bottom[j + 1]
                Q = np.array([L[u1], L[u2], newpos[b2], newpos[b1]])
                w, d, m = bilinear(Q, L[v])
                if best is None or m < best[2]:
                    best = (w, d, m, [u1, u2, merged[b2], merged[b1]], Q)
            w, d, m, ids, Q = best
            face_ties.append((v, list(zip(ids, w)), w @ Q, d, m))
    # PeB top facets: the old PVW / PeB interface quads (PeB faces made of shared nodes) + the PeB's boundary faces
    pfaces, cnt = [], defaultdict(int)
    pe_els = [c for n, c in E.values() if n == PEB]
    for c in pe_els:
        for fc in HF:
            cnt[frozenset(c[i] for i in fc)] += 1
    for c in pe_els:
        for fc in HF:
            f = [c[i] for i in fc]
            if set(f) <= S or cnt[frozenset(f)] == 1:
                pfaces.append(f)
    PC = np.array([np.mean([L[v] for v in f], 0) for f in pfaces])

    def on_peb(p):
        best = None
        for j in np.argsort(np.linalg.norm(PC - p, axis=1))[:40]:
            f = pfaces[j]; Q = np.array([L[v] for v in f])
            w, d, m = bilinear(Q, p)
            if m <= 1.0 + 1e-4 and (best is None or d < best[1]):
                best = (list(zip(f, w)), d, w @ Q)
        return best

    # 1) the PVW's new bottom nodes (merged ones keep the curve's id) and their PeB ties
    nid = max(X) + 1
    newid = {}
    for s in sorted(S):
        if s in merged:
            newid[s] = merged[s]
        else:
            newid[s] = nid; nid += 1
    ties = ''
    nt_pvw = 0
    for s in sorted(S):
        hit = on_peb(newpos[s])
        assert hit is not None and hit[1] < 0.05, (s, hit and hit[1])
        ties += lc(newid[s], [(q, w) for q, w in hit[0]]); nt_pvw += 1
    # 2) the curves' face-column nodes on the PVW's lateral face
    moved = {newid[s]: newpos[s] for s in S}
    for v, hosts, p, d, m in face_ties:
        moved[v] = p
        hosts = [(q if q not in merged else merged[q], w) for q, w in hosts]
        ties += lc(v, hosts)
    log.append(f'{len(face_ties)} curve face-column nodes tied onto the PVW lateral faces (moved up to '
               f'{max(d for *_, d, m in face_ties):.2f} mm, max |xi| {max(m for *_, m in face_ties):.2f})')
    # 3) the curves' end-ring nodes over the PeB
    npeb, dists, rejected = 0, [], []
    from curve_swell import TETS
    curve_els = {k: c for k, (n, c) in E.items() if n.startswith('canal_wrap_hex')}
    by_node = defaultdict(list)
    for k, c in curve_els.items():
        for v in c:
            by_node[v].append(c)

    def min_subtet(v, mv):
        P = lambda q: mv.get(q, L[q])   # noqa: E731
        return min(np.dot(P(c[b]) - P(c[a]), np.cross(P(c[d]) - P(c[a]), P(c[e]) - P(c[a]))) / 6
                   for c in by_node[v] for a, b, d, e in TETS)
    if mode == 'full':
        done = set(moved)
        for ring in end_rings:
            for v in ring:
                if v in done:
                    continue
                hit = on_peb(L[v])
                if hit is None or hit[1] > DTOL:
                    continue
                old = moved.get(v)
                moved[v] = hit[2]
                if min_subtet(v, moved) < 0.05:                       # the snap would crush one of its hexes: leave it free
                    rejected.append(v)
                    if old is None:
                        del moved[v]
                    else:
                        moved[v] = old
                    continue
                dists.append(hit[1])
                ties += lc(v, hit[0]); npeb += 1
        log.append(f'{len(rejected)} candidates left free (snapping would crush a hex: {rejected}); {npeb} curve end nodes snapped onto the PeB top and tied (moved {np.round(min(dists), 2) if dists else 0}'
                   f'-{np.round(max(dists), 2) if dists else 0} mm)')
    # ---- write
    Lnew = dict(L); Lnew.update(moved)
    Xg = {k: c0 + B @ p for k, p in Lnew.items()}
    for v in moved:
        if v in X:
            t = re.sub(r'(<node id="%d">)[^<]*(</node>)' % v, lambda m: m.group(1) + ','.join(f'{c:.7g}' for c in Xg[v]) + m.group(2), t, count=1)
    newnodes = ''.join(f'\t\t\t<node id="{newid[s]}">' + ','.join(f'{c:.7g}' for c in Xg[newid[s]]) + '</node>\n'
                       for s in sorted(S) if newid[s] not in X)
    a = t.index('<Nodes name="apex_contour_nodes">'); b = t.index('</Nodes>', a) + len('</Nodes>\n')
    t = t[:b] + '\t\t<Nodes name="junction_pvw_bottom">\n' + newnodes + '\t\t</Nodes>\n' + t[b:]
    remap = {s: newid[s] for s in S}
    for k, c in PL.items():
        t = t.replace(f'<elem id="{k}">' + ','.join(map(str, c)) + '</elem>',
                      f'<elem id="{k}">' + ','.join(str(remap.get(v, v)) for v in c) + '</elem>', 1)
    plf = {frozenset(c[i] for i in fc) for c in PL.values() for fc in HF}
    nsurf = defaultdict(int)

    def fix_facet(m):
        ids = [int(v) for v in m.group(2).split(',')]
        if frozenset(ids) in plf and set(ids) & S:
            nsurf[m.group(0)[:0] or 'n'] += 1
            return m.group(1) + ','.join(str(remap.get(v, v)) for v in ids) + m.group(3)
        return m.group(0)
    for sm in re.finditer(r'<Surface name="([^"]+)">.*?</Surface>', t, re.S):
        pass
    t = re.sub(r'(<(?:quad4|tri3) id="\d+">)([\d,]+)(</(?:quad4|tri3)>)', fix_facet, t)
    mid_new = [newid[s] for s in S if abs(L[s][0]) < 0.05 and newid[s] not in X]
    m = re.search(r'(<NodeSet name="BC-VW-mid">)(.*?)(</NodeSet>)', t, re.S)
    if mid_new:
        t = t[:m.start(2)] + m.group(2).rstrip() + ',\n\t\t\t' + ', '.join(map(str, mid_new)) + '\n\t\t' + t[m.end(2):]
    blk = re.search(r'\t\t<constraint name="canal_wrap_embed" type="linear constraint">.*?</constraint>\n', t, re.S).group(0)
    t = t.replace(blk, '\t\t<constraint name="junction_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n'
                  '\t\t\t<penalty>10</penalty>\n\t\t\t<maxaug>10</maxaug>\n' + ties + '\t\t</constraint>\n')
    # ---- checks
    Xc = {k: c0 + B @ p for k, p in Lnew.items()}
    for s in S:
        Xc[newid[s]] = c0 + B @ newpos[s]
    E2 = {k: (n, [remap.get(v, v) for v in c] if k in PL else c) for k, (n, c) in E.items()}
    ov = overlap(Xc, E2)
    from curve_swell import TETS
    chk = list(PL) + [k for r in sides for a, b in r[1:] for k in range(a, b + 1)]
    worst = min((np.einsum('j,j->', Xc[E2[k][1][q1]] - Xc[E2[k][1][q0]], np.cross(Xc[E2[k][1][q2]] - Xc[E2[k][1][q0]],
                 Xc[E2[k][1][q3]] - Xc[E2[k][1][q0]])) / 6, k) for k in chk for q0, q1, q2, q3 in TETS)
    log.append(f'checks: curve hexes inside PVW / PeB {len(ov)}; smallest sub-tet volume of the changed hexes {worst[0]:.3f} mm3 '
               f'(element {worst[1]}); {nt_pvw} PVW bottom nodes tied to the PeB ({len(S) - len(merged)} new + {len(merged)} curve '
               f'nodes); facets re-pointed {sum(nsurf.values())}; BC-VW-mid + {len(mid_new)}')
    print('\n'.join(log))
    assert worst[0] > 0, 'inverted sub-tet'
    name = f'{NAMES[mode]}_{case}'
    msg = (f'NOT IN SOURCE: the distal junction rebuilt (tools/build_batch278.py {mode}): the PVW flare removed (last layer '
           f'straight, own bottom nodes tied onto the PeB top), the curves joined to the straight PVW face'
           + (', the curve end rings tied onto the PeB top' if mode == 'full' else ' (curve ends NOT tied to the PeB)')
           + '; canal_wrap_embed removed.\n' + '\n'.join(log))
    print(name); print(msg)
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')


if __name__ == '__main__':
    for case in sys.argv[2:]:
        build(sys.argv[1], case)
