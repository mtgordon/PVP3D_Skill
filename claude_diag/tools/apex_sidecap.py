"""A curved cap coming out sideways from the cervix's curved end (2026-10-09, the user: "I thought I would see more of a curved cap
coming out to the side from the curved end of the cervix"; "yes, build it and show me a preview").
On each side the cervix's rounded lateral end beyond the curves (its boundary faces facing sideways, outward normal . lateral > COS,
depth beyond the curves' end + DMIN) is extruded straight out to the side (along the lateral axis) to the curves' outer lateral
extent there, in NL layers: the cap's profile is the cervix end's own curve. Its inner nodes are the cervix's surface nodes
(bonded). The curves' apical end-cap nodes (not shared with the walls) are tied to the cap's closest facet point (displacement
ties). Covered load / contact faces move to the cap's outer face. Material and damping: the AVW curve's. NOT IN SOURCE.
usage: py -3.10 apex_sidecap.py SRC NAME [--dry] [--preview PNG]"""
import re
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

from apex_close import HF, hex_jac, _closest_tri

WALLS = ('_PickedSet346', '_PickedSet64', '_PickedSet347', '_PickedSet66')
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')


def sidecap(t, c0, nn, NL=4, COS=0.5, DMIN=1.0):
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')
         if e.get('name') in WALLS + CURVES}
    H = lambda p: (p - c0) @ nn  # noqa: E731
    e3 = -nn
    e1 = np.array([1.0, 0, 0]); e1 -= (e1 @ e3) * e3; e1 /= np.linalg.norm(e1)
    cnt = Counter()
    for d in WALLS + CURVES:
        for el in E.get(d, []):
            for fc in HF:
                cnt[frozenset(el[i] for i in fc)] += 1
    curvn = [v for d in CURVES for el in E[d] for v in el]
    hcurve = min(H(X[v]) for v in curvn)
    nid = max(X) + 1
    new_nodes, hexes, outer, log = [], [], {}, []
    pos = lambda m: X[m] if m in X else new_nodes[m - nid]  # noqa: E731
    for sgn in (1, -1):
        lat = sgn * e1
        xo = max(X[v] @ lat for v in curvn if X[v] @ lat > 0 and H(X[v]) < hcurve + 3)
        faces = []
        for el in E['_PickedSet346']:
            for fc in HF:
                f = [el[i] for i in fc]
                if cnt[frozenset(f)] != 1:
                    continue
                Q = np.array([X[v] for v in f])
                nv = np.cross(Q[2] - Q[0], Q[3] - Q[1]); nv /= np.linalg.norm(nv)
                if nv @ lat > COS and max(H(X[v]) for v in f) < hcurve + DMIN:
                    faces.append(f)
        lay = {}
        for v in sorted({v for f in faces for v in f}):
            ln = xo - X[v] @ lat
            for k in range(1, NL + 1):
                lay[(v, k)] = nid + len(new_nodes)
                new_nodes.append(X[v] + lat * max(ln, 0.5) * k / NL)
        node = lambda v, k: v if k == 0 else lay[(v, k)]  # noqa: E731
        bad = 0
        for f in faces:
            for k in range(NL):
                h = [node(v, k) for v in f] + [node(v, k + 1) for v in f]
                if hex_jac(np.array([pos(m) for m in h])) <= 0:
                    h = [node(v, k) for v in f[::-1]] + [node(v, k + 1) for v in f[::-1]]
                    if hex_jac(np.array([pos(m) for m in h])) <= 0:
                        bad += 1
                hexes.append(h)
            outer[frozenset(f)] = {v: node(v, NL) for v in f}
        L = [xo - X[v] @ lat for f in faces for v in f]
        log.append(f'{"right" if sgn > 0 else "left"}: {len(faces)} cervix end faces (depth '
                   f'{min(-H(X[v]) for f in faces for v in f):.1f}..{max(-H(X[v]) for f in faces for v in f):.1f}) out to the side to '
                   f'the curves\' outer edge ({xo:.1f} mm; {min(L):.1f}..{max(L):.1f} mm long), {len(faces) * NL} hexes, {bad} inverted')
    # ties: the curves' apical end caps to the side caps
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
    wallnodes = {v for d in WALLS for el in E.get(d, []) for v in el}
    capn = set()
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
    cons, gaps = [], []
    for v in sorted(capn - wallnodes):
        dd, i, q, w = _closest_tri(X[v], TR)
        gaps.append(dd)
        cons.append((v, [(m, c) for m, c in zip(trin[i], w) if c > 1e-9]))
    # write
    blk = '\t\t<Nodes name="apex_sidecap_nodes">\n' + ''.join(
        f'\t\t\t<node id="{nid + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n' for i, p in enumerate(new_nodes)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    i = t.index('\t</Mesh>')
    t = t[:i] + '\t\t<Elements type="hex8" name="apex_sidecap">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hexes)) + '\t\t</Elements>\n' + t[i:]
    mat = re.search(r'<SolidDomain name="canal_wrap_hex_avw" mat="([^"]+)"', t).group(1)
    i = t.index('\t</MeshDomains>')
    t = t[:i] + f'\t\t<SolidDomain name="apex_sidecap" mat="{mat}" />\n' + t[i:]
    cb = re.search(r'<body_load name="wall_damping_canal_wrap_hex_avw"[^>]*>\s*<C>([^<]+)</C>', t)
    if cb:
        i = t.index('\t</Loads>')
        t = t[:i] + (f'\t\t<body_load name="wall_damping_apex_sidecap" type="mass damping" elem_set="apex_sidecap">\n'
                     f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n') + t[i:]
    cx = ('\t\t<constraint name="apex_sidecap_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
          '\t\t\t<maxaug>10</maxaug>\n' + ''.join(
              f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{c}">1</node>\n'
              + ''.join(f'\t\t\t\t<node id="{w}" bc="{c}">{-cw:.6g}</node>\n' for w, cw in ws)
              + '\t\t\t</linear_constraint>\n' for m, ws in cons for c in 'xyz') + '\t\t</constraint>\n')
    i = t.index('\t</Constraints>')
    t = t[:i] + cx + t[i:]
    moved = Counter()
    for sm in list(re.finditer(r'<Surface name="([^"]+)">', t)):
        nm = sm.group(1)
        a = t.index(f'<Surface name="{nm}">')
        b = t.index('</Surface>', a)

        def swap(mm, nm=nm):
            ids = [int(v) for v in mm.group(2).split(',')]
            mp = outer.get(frozenset(ids))
            if mp is None:
                return mm.group(0)
            moved[nm] += 1
            return f'{mm.group(1)}{",".join(str(mp[v]) for v in ids)}{mm.group(3)}'
        t = t[:a] + re.sub(r'(<(?:quad4|tri3) id="\d+">)([^<]+)(</(?:quad4|tri3)>)', swap, t[a:b]) + t[b:]
    msg = ('NOT IN SOURCE (the user: a curved cap coming out to the side from the cervix\'s curved end; tools/apex_sidecap.py): '
           + '; '.join(log) + f'; {len(new_nodes)} nodes, {len(hexes)} hexes in apex_sidecap ({mat}); {len(cons)} curve end-cap nodes '
           f'tied to the caps (rest distance median {np.median(gaps):.2f}, max {max(gaps):.2f} mm); covered faces moved to the cap\'s '
           'outer face: ' + (', '.join(f'{k} {v}' for k, v in moved.items()) or 'none'))
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
    t2, msg = sidecap(t, c0, nn)
    print(msg)
    if '--dry' not in sys.argv:
        write(name, t2, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')
