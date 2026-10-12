"""Close the tube's apex with a dome over the cervix's apical end (2026-10-09, the user: "take the very apical side of the cervix
where it rounds from the avw to the pvw and create a curved dome-like structure there that connects to the ledges [the tube's
curves here]. The dome will not be exactly spherical"). Replaces the axis extrusions of apex_close.py (v3-v5), which stopped
short of the cervix end or did not follow its shape.
How: the dome grows outward from the existing OUTER free surface of the apical end (cervix + PVW boundary faces beyond where the
curves end, H < curve end + DH; the canal-lumen faces of SlidingElastic1 excluded), so it takes the end's own contour (the round
over from the cervix front face to the PVW back face, and the lateral corners). Its inner nodes ARE the wall's surface nodes
(bonded, as the curves are); NL layers to the curves' thickness T along area-weighted nodal normals. The curve end caps (their
free apical faces) are tied to the closest point of the dome's boundary (linear constraints on displacement, penalty 10, maxaug
10; nodes shared with the walls are skipped). Every load / contact surface face covered by the dome moves to the dome's outer face
(same node order). The dome gets the curves' material (canal_wrap_hex_avw's) and their mass damping. NOT IN SOURCE.
usage (module): t2, msg = dome(t, c0, nn, T=2.3, NL=4, DH=1.0)"""
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

import numpy as np

from apex_close import HF, hex_jac, _closest_tri

WALLS = ('_PickedSet346', '_PickedSet64', '_PickedSet347', '_PickedSet66')
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
ENDWALLS = ('_PickedSet346', '_PickedSet64')


def dome(t, c0, nn, T=2.3, NL=4, DH=1.0, sides=True):
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')
         if e.get('name') in WALLS + CURVES}
    H = lambda p: (p - c0) @ nn  # noqa: E731
    cnt = Counter()
    for d in WALLS + CURVES:
        for el in E.get(d, []):
            for fc in HF:
                cnt[frozenset(el[i] for i in fc)] += 1
    lumen = set()
    for s in ('SlidingElastic1Secondary', 'SlidingElastic1Primary'):
        b = re.search(r'<Surface name="%s">(.*?)</Surface>' % s, t, re.S).group(1)
        lumen |= {frozenset(int(v) for v in x.split(',')) for x in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', b)}
    hcurve = min(H(X[v]) for d in CURVES for el in E[d] for v in el)
    # sides only (the user, 2026-10-09: "it just needs to be on the side so it can attach to the curved edges"): on each side
    # keep the faces at least as far out laterally as that side's curve end (the innermost curve node near the apex)
    xcut = {}
    for sgn in (1, -1):
        xs = [abs(X[v][0]) for d in CURVES for el in E[d] for v in el
              if sgn * X[v][0] > 1.5 and H(X[v]) < hcurve + 3.0]
        xcut[sgn] = min(xs) if (sides and xs) else 0.0
    faces = []                                              # outward-ordered quads (HF order of their hex)
    for d in ENDWALLS:
        for el in E[d]:
            for fc in HF:
                f = [el[i] for i in fc]
                k = frozenset(f)
                if cnt[k] == 1 and k not in lumen and max(H(X[v]) for v in f) < hcurve + DH:
                    xc = np.mean([X[v][0] for v in f])
                    if abs(xc) >= xcut[1 if xc > 0 else -1]:
                        faces.append(f)
    # area-weighted nodal normals over the dome faces (outward)
    nsum = defaultdict(lambda: np.zeros(3))
    for f in faces:
        P = np.array([X[v] for v in f])
        nv = np.cross(P[2] - P[0], P[3] - P[1])            # HF order: outward
        for v in f:
            nsum[v] += nv
    N = {v: s / np.linalg.norm(s) for v, s in nsum.items()}
    nid = max(X) + 1
    lay = {}
    new_nodes = []
    for v in sorted(N):
        for k in range(1, NL + 1):
            lay[(v, k)] = nid + len(new_nodes)
            new_nodes.append(X[v] + N[v] * T * k / NL)
    node = lambda v, k: v if k == 0 else lay[(v, k)]  # noqa: E731
    pos = lambda m: X[m] if m in X else new_nodes[m - nid]  # noqa: E731
    hexes, bad = [], 0
    for f in faces:
        for k in range(NL):
            h = [node(v, k) for v in f] + [node(v, k + 1) for v in f]
            if hex_jac(np.array([pos(m) for m in h])) <= 0:
                h = [node(v, k) for v in f[::-1]] + [node(v, k + 1) for v in f[::-1]]
                if hex_jac(np.array([pos(m) for m in h])) <= 0:
                    bad += 1
            hexes.append(h)
    outer = {frozenset(f): [node(v, NL) for v in f] for f in faces}
    # the dome's boundary faces (for the curve-cap ties), as triangles with node ids
    dc = Counter()
    for h in hexes:
        for fc in HF:
            dc[frozenset(h[i] for i in fc)] += 1
    tris, trin = [], []
    for h in hexes:
        for fc in HF:
            f = [h[i] for i in fc]
            if dc[frozenset(f)] == 1:
                tris += [tuple(pos(m) for m in (f[0], f[1], f[2])), tuple(pos(m) for m in (f[0], f[2], f[3]))]
                trin += [(f[0], f[1], f[2]), (f[0], f[2], f[3])]
    TR = np.array(tris)
    wallnodes = {v for d in WALLS for el in E[d] for v in el}
    capn = set()                                            # the curves' apical end caps (as apex_close.py finds them)
    for d in CURVES:
        for el in E[d]:
            hc = np.array([X[v] for v in el]).mean(0)
            for fc in HF:
                f = [el[i] for i in fc]
                if cnt[frozenset(f)] != 1:
                    continue
                Q = np.array([X[v] for v in f])
                if abs(Q.mean(0)[0]) <= 1.5 or H(Q.mean(0)) > -70:
                    continue
                nv = np.cross(Q[2] - Q[0], Q[3] - Q[1]); nv /= np.linalg.norm(nv)
                if nv @ (Q.mean(0) - hc) < 0:
                    nv = -nv
                if nv @ nn < -0.6:
                    capn |= set(f)
    ncap = len(capn)
    capn -= wallnodes
    cons, gaps = [], []
    for v in sorted(capn):
        dd, i, q, w = _closest_tri(X[v], TR)
        gaps.append(dd)
        cons.append((v, [(m, c) for m, c in zip(trin[i], w) if c > 1e-9]))
    # write: nodes, elements, domain, damping, ties, surfaces
    blk = '\t\t<Nodes name="apex_dome_nodes">\n' + ''.join(
        f'\t\t\t<node id="{nid + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n' for i, p in enumerate(new_nodes)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    eblk = '\t\t<Elements type="hex8" name="apex_dome">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hexes)) + '\t\t</Elements>\n'
    i = t.index('\t</Mesh>')
    t = t[:i] + eblk + t[i:]
    mat = re.search(r'<SolidDomain name="canal_wrap_hex_avw" mat="([^"]+)"', t).group(1)
    i = t.index('\t</MeshDomains>')
    t = t[:i] + f'\t\t<SolidDomain name="apex_dome" mat="{mat}" />\n' + t[i:]
    cb = re.search(r'<body_load name="wall_damping_canal_wrap_hex_avw"[^>]*>\s*<C>([^<]+)</C>', t)
    if cb:
        i = t.index('\t</Loads>')
        t = t[:i] + (f'\t\t<body_load name="wall_damping_apex_dome" type="mass damping" elem_set="apex_dome">\n'
                     f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n') + t[i:]
    cx = ('\t\t<constraint name="apex_dome_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
          '\t\t\t<maxaug>10</maxaug>\n' + ''.join(
              f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{c}">1</node>\n'
              + ''.join(f'\t\t\t\t<node id="{w}" bc="{c}">{-cw:.6g}</node>\n' for w, cw in ws)
              + '\t\t\t</linear_constraint>\n' for m, ws in cons for c in 'xyz') + '\t\t</constraint>\n')
    i = t.index('\t</Constraints>')
    t = t[:i] + cx + t[i:]
    moved = Counter()

    def swap(mm):
        ids = [int(v) for v in mm.group(2).split(',')]
        o = outer.get(frozenset(ids))
        if o is None:
            return mm.group(0)
        # keep the face's own node order: map each old node to its outer copy
        moved[cur[0]] += 1
        mp = {v: lay[(v, NL)] for v in ids}
        return f'{mm.group(1)}{",".join(str(mp[v]) for v in ids)}{mm.group(3)}'
    for sm in list(re.finditer(r'<Surface name="([^"]+)">', t)):
        nm = sm.group(1)
        a = t.index(f'<Surface name="{nm}">')
        b = t.index('</Surface>', a)
        cur = [nm]
        t = t[:a] + re.sub(r'(<(?:quad4|tri3) id="\d+">)([^<]+)(</(?:quad4|tri3)>)', swap, t[a:b]) + t[b:]
    where = (f'a dome-like cap on each SIDE of the cervix\'s apical end only (lateral |x| >= {xcut[1]:.1f} right, '
             f'{xcut[-1]:.1f} left), connecting the curves' if sides else 'a dome over the whole apical end')
    msg = (f'NOT IN SOURCE (the user: {where}; tools/apex_dome.py): '
           f'{len(faces)} outer free faces of the cervix / PVW beyond the curve ends (H < {hcurve + DH:.1f}) grown outward '
           f'{T:g} mm in {NL} layers along nodal normals ({len(new_nodes)} nodes, {len(hexes)} hexes in apex_dome, material {mat}, '
           f'{bad} with a non-positive Jacobian); {ncap} curve end-cap nodes, {len(cons)} of them (not shared with the walls) tied to the dome boundary (rest distance median '
           f'{np.median(gaps):.2f}, max {max(gaps):.2f} mm); covered faces moved to the dome\'s outer face: '
           + (', '.join(f'{k} {v}' for k, v in moved.items()) or 'none'))
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
    t2, msg = dome(t, c0, nn)
    print(msg)
    if '--dry' not in sys.argv:
        write(name, t2, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')
