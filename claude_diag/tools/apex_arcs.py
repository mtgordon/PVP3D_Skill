"""Close the tube's apical corners with semicircular arcs (2026-10-09, the user: "From where you stopped the curved edges, draw
semicircles up the vaginal wall (and cervix) that connect the top edge (top of the cervix) to the bottom edge (bottom of the pvw).
The radius of the semicircle will be the distance that it has to cross. Those arcs will form the outer surface of the cap to
fill in the hole." A semicircle with its ends on the two edges has a DIAMETER equal to the distance crossed: radius = half.)
At stations every ~STEP mm of depth, from the curves' apical end cap up to the cervix's tip, the top edge (the wall's lateral
extremity with the largest AP, the cervix / AVW side) and the bottom edge (smallest AP: the PVW, then the cervix's back where
the PVW has ended) are found; the arc runs from the top edge to the bottom edge, bulging out to the side, in the plane of the
lateral and AP axes. The arcs are the cap's outer surface; the cap is a shell of the curves' thickness T under them (5 x 7
nodes per station, as the curves' cap: 5 through the thickness, 7 round the arc), its first station being the curves' end cap
itself (shared nodes: the cap carries on from the curves). The arc ends are tied to the closest wall facet point (displacement
ties). New elements apex_arcs (the AVW curve's material and damping). NOT IN SOURCE.
usage: py -3.10 apex_arcs.py SRC NAME [--dry]"""
import re
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

from apex_close import HF, hex_jac, _closest_tri
from apex_seal import _grid

CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
WALLS = ('_PickedSet346', '_PickedSet64', '_PickedSet347')


def arcs(t, c0, nn, STEP=0.7, T=2.3, MINSPAN=1.5, const_r=False):
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')
         if e.get('name') in CURVES + WALLS}
    H = lambda p: (p - c0) @ nn  # noqa: E731
    e3 = -nn
    e1 = np.array([1.0, 0, 0]); e1 -= (e1 @ e3) * e3; e1 /= np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    B = np.array([e1, e2, e3]).T
    L = lambda p: (p - c0) @ B  # noqa: E731
    G = lambda q: c0 + np.asarray(q) @ B.T  # noqa: E731
    cnt = Counter()
    for d in CURVES:
        for el in E[d]:
            for fc in HF:
                cnt[frozenset(el[i] for i in fc)] += 1
    wc = Counter()
    for d in WALLS:
        for el in E[d]:
            for fc in HF:
                wc[frozenset(el[i] for i in fc)] += 1
    tris, trin = [], []
    for d in WALLS:
        for el in E[d]:
            for fc in HF:
                f = [el[i] for i in fc]
                if wc[frozenset(f)] == 1:
                    tris += [(X[f[0]], X[f[1]], X[f[2]]), (X[f[0]], X[f[2]], X[f[3]])]
                    trin += [(f[0], f[1], f[2]), (f[0], f[2], f[3])]
    TR = np.array(tris)
    wallv = {v for d in WALLS for el in E[d] for v in el}
    nid = max(X) + 1
    new_nodes, hexes, cons, log = [], [], [], []
    pos = lambda m: X[m] if m in X else new_nodes[m - nid]  # noqa: E731

    def add(q):
        new_nodes.append(G(q))
        return nid + len(new_nodes) - 1
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
        Gr = np.array(_grid([f for d, f in caps]))
        if Gr.shape[0] > Gr.shape[1]:
            Gr = Gr.T
        assert Gr.shape == (5, 7), Gr.shape
        if 'canal_wrap_hex_avw' not in owner[Gr[2, 0]]:
            Gr = Gr[:, ::-1]                                # column 0 = top (cervix / AVW) side
        cen = np.array([L(X[v]) for v in Gr.ravel()]).mean(0)
        rad = lambda r: np.mean([np.hypot(L(X[v])[1] - cen[1], L(X[v])[0] - cen[0]) for v in Gr[r]])  # noqa: E731
        if rad(0) > rad(4):
            Gr = Gr[::-1, :]                                # row 0 = inner (smaller radius)
        d0 = max(L(X[v])[2] for v in Gr.ravel())            # the cap's deepest point
        # the true edges (2026-10-09): where the canal's lumen surfaces end at this corner. Top edge = the boundary of the
        # cervix / AVW lumen surface (SlidingElastic1Secondary), bottom edge = the PVW's (SlidingElastic1Primary); their lateral
        # runs above the curves' end, interpolated by depth
        def edge(sname):
            b = re.search(r'<Surface name="%s">(.*?)</Surface>' % sname, t, re.S).group(1)
            F = [[int(v) for v in x.split(',')] for x in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', b)]
            ec = Counter()
            for f in F:
                for k in range(len(f)):
                    ec[tuple(sorted((f[k], f[(k + 1) % len(f)])))] += 1
            nodes = {v for e, c in ec.items() if c == 1 for v in e if sgn * X[v][0] > 0}
            P = np.array(sorted([L(X[v]) for v in nodes if L(X[v])[2] > d0 - 1.0], key=lambda p: p[2]))
            xc = sgn * max(sgn * P[:, 0])
            P = P[np.abs(P[:, 0] - xc) < 1.5]               # the lateral run only (not the run across the top)
            return P
        Pt, Pb = edge('SlidingElastic1Secondary'), edge('SlidingElastic1Primary')
        dend = min(Pt[:, 2].max(), Pb[:, 2].max())
        nst = max(2, int(round((dend - d0) / STEP)) + 1)
        stations = []
        for dd in np.linspace(d0 + (dend - d0) / nst, dend, nst):
            top = np.array([np.interp(dd, Pt[:, 2], Pt[:, 0]), np.interp(dd, Pt[:, 2], Pt[:, 1]), dd])
            bot = np.array([np.interp(dd, Pb[:, 2], Pb[:, 0]), np.interp(dd, Pb[:, 2], Pb[:, 1]), dd])
            stations.append((dd, top, bot))
        # the user (2026-10-09): "keep going even though the pvw is not there on the other side": beyond the PVW's end the arcs
        # carry on up the cervix to its tip with the same ends (same lateral / AP positions) and the same radius
        dtip = max(L(X[v])[2] for el in E['_PickedSet346'] for v in el if sgn * X[v][0] > 0)
        _, top_l, bot_l = stations[-1]
        dd = dend + STEP
        while dd < dtip - 0.2:
            stations.append((dd, np.array([top_l[0], top_l[1], dd]), np.array([bot_l[0], bot_l[1], dd])))
            dd += STEP
        stations.append((dtip, np.array([top_l[0], top_l[1], dtip]), np.array([bot_l[0], bot_l[1], dtip])))
        rings = [Gr]
        Rc = max(s[1][1] - s[2][1] for s in stations) / 2   # constant radius (the user, 2026-10-09): half the widest span
        for dd, top, bot in stations:
            # the arc in the (lateral, AP) plane at depth dd: ends at the top / bottom edges, bulging out to the side
            P_t, P_b = np.array(top[:2]), np.array(bot[:2])
            mid = (P_t + P_b) / 2
            u = P_t - P_b; half = np.linalg.norm(u) / 2; u /= np.linalg.norm(u)
            lat = np.array([sgn, 0.0]); lat -= (lat @ u) * u; lat /= np.linalg.norm(lat)
            if const_r:
                # the smallest arc of radius Rc between the two edges: centre inward of the chord, arc on the outer side
                R0 = max(Rc, half)
                cen = mid - lat * np.sqrt(max(R0 ** 2 - half ** 2, 0.0))
                a_t = np.arctan2((P_t - cen) @ lat, (P_t - cen) @ u)
                a_b = np.arctan2((P_b - cen) @ lat, (P_b - cen) @ u)
                angs = [a_t + (a_b - a_t) * b / 6 for b in range(7)]
                thick = min(T, 0.8 * R0)
                Rk = np.zeros_like(Gr)
                for a in range(5):
                    r = R0 - thick + thick * a / 4
                    for b in range(7):
                        p2 = cen + r * (np.cos(angs[b]) * u + np.sin(angs[b]) * lat)
                        Rk[a, b] = add([p2[0], p2[1], dd])
                if len(rings) == 1:
                    # blend the first arc half-way from the curves' end cap (its shape moved to this depth) into the arc
                    for a in range(5):
                        for b in range(7):
                            q_arc = L(new_nodes[Rk[a, b] - nid]); q_cap = L(X[Gr[a, b]])
                            q = 0.5 * q_arc + 0.5 * np.array([q_cap[0], q_cap[1], dd])
                            new_nodes[Rk[a, b] - nid] = G(q)
                rings.append(Rk)
                continue
            r_out = half
            thick = min(T, 0.8 * r_out)
            Rk = np.zeros_like(Gr)
            for a in range(5):
                r = r_out - thick + thick * a / 4               # row 0 inner, row 4 = the arc (outer surface)
                for b in range(7):
                    th = np.pi * b / 6                           # 0 = top edge, pi = bottom edge
                    p2 = mid + r * (np.cos(th) * u + np.sin(th) * lat)
                    Rk[a, b] = add([p2[0], p2[1], dd])
            rings.append(Rk)
        bad = 0
        for k in range(1, len(rings)):
            for a in range(4):
                for b in range(6):
                    lo = [rings[k - 1][a, b], rings[k - 1][a, b + 1], rings[k - 1][a + 1, b + 1], rings[k - 1][a + 1, b]]
                    hi = [rings[k][a, b], rings[k][a, b + 1], rings[k][a + 1, b + 1], rings[k][a + 1, b]]
                    h = lo + hi
                    if hex_jac(np.array([pos(m) for m in h])) <= 0:
                        h = lo[::-1] + hi[::-1]
                        if hex_jac(np.array([pos(m) for m in h])) <= 0:
                            bad += 1
                    hexes.append(h)
            for a in range(5):
                for b in (0, 6):
                    m = rings[k][a, b]
                    dd_, i, q, w = _closest_tri(pos(m), TR)
                    if dd_ < 3.0:
                        cons.append((m, [(x, c) for x, c in zip(trin[i], w) if c > 1e-9]))
        spans = [s[1][1] - s[2][1] for s in stations]
        log.append(f'{"right" if sgn > 0 else "left"}: {len(stations)} arcs from depth {d0:.1f} to {stations[-1][0]:.1f} '
                   f'(span {max(spans):.1f} -> {min(spans):.1f} mm, ' + (f'constant radius {Rc:.1f} mm, the smallest arc' if const_r else 'radius = half') + f'), {(len(rings) - 1) * 24} hexes, {bad} inverted')
    blk = '\t\t<Nodes name="apex_arcs_nodes">\n' + ''.join(
        f'\t\t\t<node id="{nid + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n' for i, p in enumerate(new_nodes)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    i = t.index('\t</Mesh>')
    t = t[:i] + '\t\t<Elements type="hex8" name="apex_arcs">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hexes)) + '\t\t</Elements>\n' + t[i:]
    mat = re.search(r'<SolidDomain name="canal_wrap_hex_avw" mat="([^"]+)"', t).group(1)
    i = t.index('\t</MeshDomains>')
    t = t[:i] + f'\t\t<SolidDomain name="apex_arcs" mat="{mat}" />\n' + t[i:]
    cb = re.search(r'<body_load name="wall_damping_canal_wrap_hex_avw"[^>]*>\s*<C>([^<]+)</C>', t)
    if cb:
        i = t.index('\t</Loads>')
        t = t[:i] + (f'\t\t<body_load name="wall_damping_apex_arcs" type="mass damping" elem_set="apex_arcs">\n'
                     f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n') + t[i:]
    cx = ('\t\t<constraint name="apex_arcs_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
          '\t\t\t<maxaug>10</maxaug>\n' + ''.join(
              f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{c}">1</node>\n'
              + ''.join(f'\t\t\t\t<node id="{w}" bc="{c}">{-cw:.6g}</node>\n' for w, cw in ws)
              + '\t\t\t</linear_constraint>\n' for m, ws in cons for c in 'xyz') + '\t\t</constraint>\n')
    i = t.index('\t</Constraints>')
    t = t[:i] + cx + t[i:]
    msg = ('NOT IN SOURCE (the user: semicircular arcs from the cervix\'s top edge to the PVW\'s bottom edge, up from the curves\' end, '
           f'as the cap\'s outer surface; tools/apex_arcs.py): ' + '; '.join(log) + f'; {len(new_nodes)} nodes, {len(hexes)} hexes in '
           f'apex_arcs ({mat}, shell {T:g} mm under the arcs, first station = the curves\' end cap), {len(cons)} arc-end ties')
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
    t2, msg = arcs(t, c0, nn, const_r='--const' in sys.argv)
    print(msg)
    if '--dry' not in sys.argv:
        write(name, t2, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')
