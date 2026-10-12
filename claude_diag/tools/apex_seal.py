"""Seal the tube's apical corners (2026-10-09, the user, after the caps / dome: "Forget about the caps ... the gap showing between
the cervix and the yellow pvw ... is the shape that needs to be covered ... extend the curved ledges that form the sides and then
cap them with a shape that covers the back of the cervix"; plan agreed: "sure").
At each lateral corner the curves (canal_wrap_hex_avw = the U's front half on the cervix / AVW, canal_wrap_hex_pvw = its back
half on the PVW) stop at depth ~86-87 mm while the cervix runs to ~93 and the PVW to ~90, so the canal's side is open above them.
1. Extension: the curves' apical end cap (a 5 x 7 node grid: 5 through the thickness, 7 round the U) is swept up to the PVW's end
   (DEND, default: the deepest PVW node on that side) in K rings. Every node's depth goes linearly to DEND; its lateral / AP
   position follows the walls' lateral edges (straight-line fits over the cap-to-DEND depth range): the front arm follows the
   cervix's edge, the back arm the PVW's, the nodes between blend by their place round the U. The rings' arm tips are tied
   to the closest cervix / PVW facet point (displacement ties: the rest offset is kept).
2. Lid: the last two rings' INNER layer (the U's inside, against the walls) closes with a plate across the U's opening, from the
   front tip to the back tip: hexes between the inner U nodes and a chord, wedges (penta6) at the two tips. The chord nodes
   are tied to the closest cervix / PVW facet point.
New elements: <curve>_seal (the curve halves' materials / damping), apex_seal_lid (penta6 + hex8, the AVW curve's material).
NOT IN SOURCE.
usage (module): t2, msg = seal(t, c0, nn, K=4)"""
import re
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

from apex_close import HF, hex_jac, _closest_tri

CURVES = {'canal_wrap_hex_avw': 'front', 'canal_wrap_hex_pvw': 'back'}
WALLS = ('_PickedSet346', '_PickedSet64', '_PickedSet347')


def _grid(faces):
    """recover the structured node grid of a quad patch: list of rows"""
    adj = {}
    for f in faces:
        for k in range(4):
            a, b = f[k], f[(k + 1) % 4]
            adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    deg = Counter(v for f in faces for v in f)
    c = min(v for v in deg if deg[v] == 1)
    n1, n2 = sorted(adj[c])
    def walk(a, b):
        line = [a, b]
        while True:
            # the next node along the line: neighbour of b, not in line, sharing no face with line[-2] except through b
            cand = [x for x in adj[b] if x not in line and not any({line[-2], x} <= set(f) for f in faces)]
            if not cand:
                return line
            line.append(cand[0]); b = cand[0]
    col = walk(c, n2)
    rows = []
    prev_row = None
    for i, s in enumerate(col):
        if i == 0:
            row = walk(c, n1)
        else:
            row = [s]
            for j in range(1, len(rows[0])):
                # the node sharing a face with rows[-1][j-1], rows[-1][j] and row[j-1]
                f = [f for f in faces if {rows[-1][j - 1], rows[-1][j], row[j - 1]} <= set(f)][0]
                row.append([v for v in f if v not in (rows[-1][j - 1], rows[-1][j], row[j - 1])][0])
        rows.append(row)
    return rows


def seal(t, c0, nn, K=4, maxtie=3.0):
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')
         if e.get('name') in tuple(CURVES) + WALLS}
    H = lambda p: (p - c0) @ nn  # noqa: E731
    e3 = -nn
    e1 = np.array([1.0, 0, 0]); e1 -= (e1 @ e3) * e3; e1 /= np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    B = np.array([e1, e2, e3]).T
    L = lambda p: (p - c0) @ B  # noqa: E731  (lateral, AP, depth)
    G = lambda q: c0 + q @ B.T  # noqa: E731
    cnt = Counter()
    for d in CURVES:
        for el in E[d]:
            for fc in HF:
                cnt[frozenset(el[i] for i in fc)] += 1
    # wall boundary triangles for ties
    wc = Counter()
    for d in WALLS:
        for el in E[d]:
            for fc in HF:
                wc[frozenset(el[i] for i in fc)] += 1
    TRI = {'front': ([], []), 'back': ([], [])}
    for d, side in (('_PickedSet346', 'front'), ('_PickedSet347', 'front'), ('_PickedSet64', 'back')):
        for el in E[d]:
            for fc in HF:
                f = [el[i] for i in fc]
                if wc[frozenset(f)] == 1:
                    TRI[side][0].extend([(X[f[0]], X[f[1]], X[f[2]]), (X[f[0]], X[f[2]], X[f[3]])])
                    TRI[side][1].extend([(f[0], f[1], f[2]), (f[0], f[2], f[3])])
    TR = {k: (np.array(v[0]), v[1]) for k, v in TRI.items()}
    nid = max(X) + 1
    new_nodes, cons, log = [], [], []
    hexes = {d + '_seal': [] for d in CURVES}
    lid_hex, lid_pen = [], []
    pos = lambda m: X[m] if m in X else new_nodes[m - nid]  # noqa: E731

    def add(p):
        new_nodes.append(p)
        return nid + len(new_nodes) - 1

    def tie(m, side):
        dd, i, q, w = _closest_tri(pos(m), TR[side][0])
        if dd <= maxtie:
            cons.append((m, [(a, c) for a, c in zip(TR[side][1][i], w) if c > 1e-9]))
        return dd
    for sgn in (1, -1):
        caps = []
        for d in CURVES:
            for el in E[d]:
                hc = np.array([X[v] for v in el]).mean(0)
                for fc in HF:
                    f = [el[i] for i in fc]
                    if cnt[frozenset(f)] != 1:
                        continue
                    Q = np.array([X[v] for v in f])
                    if sgn * Q.mean(0)[0] <= 1.5 or H(Q.mean(0)) > -70:
                        continue
                    nv = np.cross(Q[2] - Q[0], Q[3] - Q[1]); nv /= np.linalg.norm(nv)
                    if nv @ (Q.mean(0) - hc) < 0:
                        nv = -nv
                    if nv @ nn < -0.6:
                        caps.append((d, f))
        owner = {}
        for d, f in caps:
            for v in f:
                owner.setdefault(v, set()).add(d)
        rows = _grid([f for d, f in caps])
        Gr = np.array(rows)                                 # shape (a, b)
        if Gr.shape[0] > Gr.shape[1]:
            Gr = Gr.T                                       # rows = through (5), columns = round the U (7)
        assert Gr.shape == (5, 7), Gr.shape
        # orient: column 0 = front (cervix / AVW half), row 0 = the inner layer (closest to the walls)
        if 'canal_wrap_hex_avw' not in owner[Gr[2, 0]]:
            Gr = Gr[:, ::-1]
        wallpts = np.array([X[v] for d in WALLS for el in E[d] for v in el if sgn * X[v][0] > 15])
        dist = lambda r: np.mean([np.min(np.linalg.norm(wallpts - X[v], axis=1)) for v in Gr[r]])  # noqa: E731
        if dist(4) < dist(0):
            Gr = Gr[::-1, :]
        # edge fits (lateral, AP) vs depth for the cervix (front) and the PVW (back) on this side
        dcap = np.array([[L(X[v])[2] for v in row] for row in Gr])
        side_pv = [v for el in E['_PickedSet64'] for v in el if sgn * X[v][0] > 0]
        dend = max(L(X[v])[2] for v in side_pv)
        fits = {}
        for nm, S in (('front', [v for el in E['_PickedSet346'] for v in el if sgn * X[v][0] > 0]), ('back', side_pv)):
            pts = []
            for dd in np.arange(dcap.min() - 2, dend + 0.01, 0.5):
                band = [L(X[v]) for v in S if abs(L(X[v])[2] - dd) < 0.5]
                if band:
                    xm = max(abs(p[0]) for p in band)
                    pts.append(np.mean([p for p in band if abs(p[0]) > xm - 1.0], axis=0))
            P = np.array(pts)
            fits[nm] = (np.polyfit(P[:, 2], P[:, 0], 1), np.polyfit(P[:, 2], P[:, 1], 1))
        ev = lambda nm, d: np.array([np.polyval(fits[nm][0], d), np.polyval(fits[nm][1], d)])  # noqa: E731
        rings = [Gr]
        for k in range(1, K + 1):
            Rk = np.zeros_like(Gr)
            for a in range(5):
                for b in range(7):
                    v = Gr[a, b]
                    q = L(X[v])
                    d0 = q[2]
                    dk = d0 + (dend - d0) * k / K
                    s = b / 6.0
                    shift = (1 - s) * (ev('front', dk) - ev('front', d0)) + s * (ev('back', dk) - ev('back', d0))
                    Rk[a, b] = add(G(np.array([q[0] + shift[0], q[1] + shift[1], dk])))
            rings.append(Rk)
        # hexes between consecutive rings, from the cap faces' pattern
        bad = 0
        for k in range(1, K + 1):
            for a in range(4):
                for b in range(6):
                    half = 'canal_wrap_hex_avw' if b < 3 else 'canal_wrap_hex_pvw'
                    lo = [rings[k - 1][a, b], rings[k - 1][a, b + 1], rings[k - 1][a + 1, b + 1], rings[k - 1][a + 1, b]]
                    hi = [rings[k][a, b], rings[k][a, b + 1], rings[k][a + 1, b + 1], rings[k][a + 1, b]]
                    h = lo + hi
                    if hex_jac(np.array([pos(m) for m in h])) <= 0:
                        h = lo[::-1] + hi[::-1]
                        if hex_jac(np.array([pos(m) for m in h])) <= 0:
                            bad += 1
                    hexes[half + '_seal'].append(h)
            for a in range(5):                              # arm tips follow the walls
                tie(rings[k][a, 0], 'front')
                tie(rings[k][a, 6], 'back')
        # lid between the last two rings' inner layers: hexes from the inner U to a chord, wedges at the tips
        lidbad = 0
        chord = {}
        for k in (K - 1, K):
            I = [rings[k][0, b] for b in range(7)]
            a0, a6 = pos(I[0]), pos(I[6])
            chord[k] = [I[0]] + [add(a0 + (a6 - a0) * b / 6.0) for b in range(1, 6)] + [I[6]]
        for b in range(6):
            lo_I, hi_I = rings[K - 1][0, :], rings[K][0, :]
            lo_C, hi_C = chord[K - 1], chord[K]
            if b == 0 or b == 5:
                # wedge: triangle (I_b, I_b+1, C_inner) on both layers
                ci = 1 if b == 0 else 5
                tri_lo = [lo_I[b], lo_I[b + 1], lo_C[ci]]
                tri_hi = [hi_I[b], hi_I[b + 1], hi_C[ci]]
                P = np.array([pos(m) for m in tri_lo + tri_hi])
                if np.cross(P[1] - P[0], P[2] - P[0]) @ (P[3:].mean(0) - P[:3].mean(0)) < 0:
                    tri_lo, tri_hi = tri_lo[::-1], tri_hi[::-1]      # penta6: the bottom triangle's normal toward the top
                lid_pen.append(tri_lo + tri_hi)
            else:
                lo = [lo_I[b], lo_I[b + 1], lo_C[b + 1], lo_C[b]]
                hi = [hi_I[b], hi_I[b + 1], hi_C[b + 1], hi_C[b]]
                h = lo + hi
                if hex_jac(np.array([pos(m) for m in h])) <= 0:
                    h = lo[::-1] + hi[::-1]
                    if hex_jac(np.array([pos(m) for m in h])) <= 0:
                        lidbad += 1
                lid_hex.append(h)
        tied = 0
        for k in (K - 1, K):
            for b in range(1, 6):
                side = 'front' if b < 3 else 'back'
                if tie(chord[k][b], side) <= maxtie:
                    tied += 1
        log.append(f'{"right" if sgn > 0 else "left"}: cap 5 x 7 at depth {dcap.min():.1f}..{dcap.max():.1f} swept to the PVW end '
                   f'(depth {dend:.1f}) in {K} rings along the cervix / PVW edges ({4 * 6 * K} hexes, {bad} inverted); lid '
                   f'{4} hexes ({lidbad} inverted) + 2 wedges, {tied} of 10 chord nodes tied')
    # write
    blk = '\t\t<Nodes name="apex_seal_nodes">\n' + ''.join(
        f'\t\t\t<node id="{nid + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n' for i, p in enumerate(new_nodes)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    eblk, dom, damp = '', '', ''
    mats = {d.get('name'): d.get('mat') for d in ET.fromstring(t).find('MeshDomains')}
    for dn, hx in list(hexes.items()) + [('apex_seal_lid', lid_hex)]:
        eblk += f'\t\t<Elements type="hex8" name="{dn}">\n' + ''.join(
            f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hx)) + '\t\t</Elements>\n'
        eid += len(hx)
        base = dn[:-len('_seal')] if dn.endswith('_seal') else 'canal_wrap_hex_avw'
        dom += f'\t\t<SolidDomain name="{dn}" mat="{mats[base]}" />\n'
        cb = re.search(r'<body_load name="wall_damping_%s"[^>]*>\s*<C>([^<]+)</C>' % re.escape(base), t)
        if cb:
            damp += (f'\t\t<body_load name="wall_damping_{dn}" type="mass damping" elem_set="{dn}">\n'
                     f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n')
    eblk += '\t\t<Elements type="penta6" name="apex_seal_lid_wedge">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, w))}</elem>\n' for j, w in enumerate(lid_pen)) + '\t\t</Elements>\n'
    dom += f'\t\t<SolidDomain name="apex_seal_lid_wedge" mat="{mats["canal_wrap_hex_avw"]}" />\n'
    cb = re.search(r'<body_load name="wall_damping_canal_wrap_hex_avw"[^>]*>\s*<C>([^<]+)</C>', t)
    if cb:
        damp += (f'\t\t<body_load name="wall_damping_apex_seal_lid_wedge" type="mass damping" elem_set="apex_seal_lid_wedge">\n'
                 f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n')
    i = t.index('\t</Mesh>')
    t = t[:i] + eblk + t[i:]
    i = t.index('\t</MeshDomains>')
    t = t[:i] + dom + t[i:]
    if damp:
        i = t.index('\t</Loads>')
        t = t[:i] + damp + t[i:]
    cx = ('\t\t<constraint name="apex_seal_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
          '\t\t\t<maxaug>10</maxaug>\n' + ''.join(
              f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{c}">1</node>\n'
              + ''.join(f'\t\t\t\t<node id="{w}" bc="{c}">{-cw:.6g}</node>\n' for w, cw in ws)
              + '\t\t\t</linear_constraint>\n' for m, ws in cons for c in 'xyz') + '\t\t</constraint>\n')
    i = t.index('\t</Constraints>')
    t = t[:i] + cx + t[i:]
    msg = ('NOT IN SOURCE (the user: extend the curves along the cervix / PVW edges to the PVW end and lid the U over the back of '
           f'the cervix; tools/apex_seal.py): {len(new_nodes)} nodes, {sum(len(v) for v in hexes.values())} sweep hexes, lid '
           f'{len(lid_hex)} hexes + {len(lid_pen)} wedges, {len(cons)} node ties (apex_seal_ties): ' + '; '.join(log))
    return t, msg


if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from variants7 import JOBS
    from build_batch162 import write
    from build_batch202 import plane
    RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
    src, name = sys.argv[1], sys.argv[2]
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    c0, nrm = plane()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    cxn = [int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet346' for el in e for v in el.text.split(',')]
    nn = nrm if np.dot(np.mean([X[v] for v in cxn], 0) - c0, nrm) < 0 else -nrm
    t2, msg = seal(t, c0, nn)
    print(msg)
    if '--dry' not in sys.argv:
        write(name, t2, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')
