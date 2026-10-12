"""Preview (no model written), 2026-10-09 after the user's "there's still a hole at the end of the U. Also, it doesn't contour
to the cervix quite as nicely as I'd like": the curves' U carried on up the canal, but every layer wraps the walls' lateral
profile AT ITS OWN DEPTH. Sections of the cervix + PVW outer surfaces give, per depth, the lateral face's position (x_w) and AP
extent (y_lo..y_hi, the 1 mm strip next to the most lateral point). Each cap node keeps its relative AP position s across that
extent and its lateral offset h from the face; the bulge shrinks with the extent (not below x0.45, so the shell keeps ~1 mm) as
the cervix rounds off to its tip, so the U narrows and closes onto the tip; over the last CLOSE layers the U's inner surface comes down onto the walls' face, so
the end is solid (no opening in the U's mouth). Output: JSON for apex_view (part 'proposal').
usage: py -3.10 apex_contour.py SRC NAME [--dry]   |   py -3.10 apex_contour.py --preview RUN OUT.json [ZEND]"""
import json
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import re  # noqa: E402
from apex_close import HF, hex_jac, _closest_tri  # noqa: E402
from apex_seal import _grid  # noqa: E402
from build_batch202 import plane  # noqa: E402
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
WALLS = ('_PickedSet346', '_PickedSet64')
CLOSE = 4          # layers over which the U's mouth closes onto the cervix at the end


def section(T, z):
    pts = []
    for tr in T:
        h = tr[:, 2] - z
        for i in range(3):
            a, b, ha, hb = tr[i], tr[(i + 1) % 3], h[i], h[(i + 1) % 3]
            if ha * hb < 0:
                pts.append(a + (b - a) * ha / (ha - hb))
    return np.array(pts)


def sweep(t, zend=None):
    """per side: (Gr node ids 5 x 7, row 0 = inner; rings (n+1, 5, 7, 3) in the canal frame, ring 0 = the cap); c0, B, X"""
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')}
    c0, nrm = plane()
    cxn = [v for el in E['_PickedSet346'] for v in el]
    nn = nrm if np.dot(np.mean([X[v] for v in cxn], 0) - c0, nrm) < 0 else -nrm
    e3 = -nn
    e1 = np.array([1.0, 0, 0]); e1 -= (e1 @ e3) * e3; e1 /= np.linalg.norm(e1)
    B = np.array([e1, np.cross(e3, e1), e3]).T
    L = {k: (p - c0) @ B for k, p in X.items()}
    cnt = Counter(frozenset(el[i] for i in fc) for d in CURVES for el in E[d] for fc in HF)
    wc = Counter(frozenset(el[i] for i in fc) for d in WALLS for el in E[d] for fc in HF)
    T = []
    for d in WALLS:
        for el in E[d]:
            for fc in HF:
                f = [el[i] for i in fc]
                if wc[frozenset(f)] == 1:
                    T += [(L[f[0]], L[f[1]], L[f[2]]), (L[f[0]], L[f[2]], L[f[3]])]
    T = np.array(T)
    T = T[T[:, :, 2].max(1) > 80]
    sides = []
    for sgn in (1, -1):
        caps = []
        for d in CURVES:
            for el in E[d]:
                hc = np.mean([L[v] for v in el], 0)
                for fc in HF:
                    f = [el[i] for i in fc]
                    if cnt[frozenset(f)] != 1:
                        continue
                    Q = np.array([L[v] for v in f])
                    if sgn * Q.mean(0)[0] <= 1.5 or Q.mean(0)[2] < 70:
                        continue
                    nv = np.cross(Q[2] - Q[0], Q[3] - Q[1]); nv /= np.linalg.norm(nv)
                    if nv @ (Q.mean(0) - hc) < 0:
                        nv = -nv
                    if nv[2] > 0.6:
                        caps.append(f)
        Gr = np.array(_grid(caps))
        if Gr.shape[0] > Gr.shape[1]:
            Gr = Gr.T
        P0 = np.array([[L[v] for v in r] for r in Gr])
        if np.abs(P0[0, 3, :2] - P0[0, 0, :2]).sum() > np.abs(P0[4, 3, :2] - P0[4, 0, :2]).sum():
            P0 = P0[::-1]; Gr = Gr[::-1]                                   # row 0 = the U's inner surface (smaller arc)
        # the walls' lateral profile by depth
        zs = np.arange(P0[:, :, 2].min() - 0.5, 93.6, 0.25)
        prof = []
        for z in zs:
            s = section(T, z)
            s = s[sgn * s[:, 0] > 10] if len(s) else s
            if len(s) < 2:
                break
            xm = (sgn * s[:, 0]).max(); t = s[sgn * s[:, 0] > xm - 1.0]
            prof.append((z, sgn * xm, t[:, 1].min(), t[:, 1].max()))
        prof = np.array(prof)
        Z, XW, YL, YH = prof.T
        f = lambda z: (np.interp(z, Z, XW), np.interp(z, Z, YL), np.interp(z, Z, YH))  # noqa: E731
        ze = zend or Z[np.where(YH - YL > 2.5)[0].max()]
        xw0, yl0, yh0 = f(P0[:, :, 2])
        S = (P0[:, :, 1] - yl0) / (yh0 - yl0)
        Hh = P0[:, :, 0] - xw0
        span0 = (yh0 - yl0)
        n = max(2, int(np.ceil((ze - P0[:, :, 2].min()) / 0.6)))
        R = np.zeros((n + 1, 5, 7, 3))
        for k in range(n + 1):
            for i in range(5):
                for j in range(7):
                    z = P0[i, j, 2] + (ze - P0[i, j, 2]) * k / n
                    xw, yl, yh = f(z)
                    g = np.clip((yh - yl) / span0[i, j], 0.45, 1.0)
                    # closing the U's mouth: over the last CLOSE layers the inner surface comes down onto the walls' face
                    c = np.clip((k - (n - CLOSE)) / CLOSE, 0, 1)
                    hh = Hh[i, j] - c * (1 - i / 4) * Hh[0, j]
                    yc = (yl + yh) / 2
                    p = np.array([xw + hh * g, yc + (yl + S[i, j] * (yh - yl) - yc) , z])
                    if k == 0:
                        p = P0[i, j]
                    R[k, i, j] = p
        msg = (f'{"right" if sgn > 0 else "left"}: cap depth {P0[:, :, 2].min():.1f}..{P0[:, :, 2].max():.1f} -> end {ze:.2f}, '
               f'AP extent {span0.mean():.1f} -> {f(ze)[2] - f(ze)[1]:.1f} mm, {n} layers')
        sides.append((Gr, R, msg))
    return sides, c0, B, X


def preview(run, out, zend=None):
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    sides, c0, B, X = sweep(t, zend)
    nodes, quads = [], []
    for Gr, R, msg in sides:
        n = len(R) - 1
        base = len(nodes)
        nodes += list(R.reshape(-1, 3))
        idx = lambda k, i, j: base + (k * 5 + i) * 7 + j  # noqa: E731
        for k in range(n):
            for i in (0, 4):
                for j in range(6):
                    quads.append([idx(k, i, j), idx(k, i, j + 1), idx(k + 1, i, j + 1), idx(k + 1, i, j)])
            for j in (0, 6):
                for i in range(4):
                    quads.append([idx(k, i, j), idx(k, i + 1, j), idx(k + 1, i + 1, j), idx(k + 1, i, j)])
        for i in range(4):
            for j in range(6):
                quads.append([idx(n, i, j), idx(n, i, j + 1), idx(n, i + 1, j + 1), idx(n, i + 1, j)])
        print(msg)
    json.dump({'nodes': [[round(float(c), 3) for c in p] for p in nodes], 'parts': {'proposal': quads}, 'edges': {}},
              open(out, 'w'))


def build(t, zend=None, TIE=1.0):
    """the cap as a model part: hexes between consecutive rings (ring 0 = the curves' own end-cap nodes, shared), domain
    apex_contour (the AVW curve's material and damping); the legs (columns 0 / 6, every row) and every inner-row node within
    TIE mm of a cervix / PVW facet are tied to their closest facet point (linear constraints, as apex_arcs)"""
    sides, c0, B, X = sweep(t, zend)
    G = lambda q: c0 + np.asarray(q) @ B.T  # noqa: E731
    mesh = ET.fromstring(t).find('Mesh')
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements') if e.get('name') in WALLS}
    wc = Counter(frozenset(el[i] for i in fc) for d in WALLS for el in E[d] for fc in HF)
    tris, trin = [], []
    for d in WALLS:
        for el in E[d]:
            for fc in HF:
                f = [el[i] for i in fc]
                if wc[frozenset(f)] == 1 and max(abs(X[v][0]) for v in f) > 10:
                    tris += [(X[f[0]], X[f[1]], X[f[2]]), (X[f[0]], X[f[2]], X[f[3]])]
                    trin += [(f[0], f[1], f[2]), (f[0], f[2], f[3])]
    TR = np.array(tris)
    nid = max(X) + 1
    new_nodes, hexes, cons, log = [], [], [], []
    pos = lambda m: X[m] if m in X else new_nodes[m - nid]  # noqa: E731
    for Gr, R, msg in sides:
        rings = [Gr]
        for k in range(1, len(R)):
            Rk = np.zeros_like(Gr)
            for i in range(5):
                for j in range(7):
                    new_nodes.append(G(R[k, i, j])); Rk[i, j] = nid + len(new_nodes) - 1
            rings.append(Rk)
        bad, nt = 0, 0
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
                for b in range(7):
                    if not (b in (0, 6) or a == 0):
                        continue
                    m = rings[k][a, b]
                    dd, i, q, w = _closest_tri(pos(m), TR)
                    if (b in (0, 6) and dd < 3.0) or dd < TIE:
                        cons.append((m, [(x, c) for x, c in zip(trin[i], w) if c > 1e-9])); nt += 1
        log.append(f'{msg}, {(len(rings) - 1) * 24} hexes, {bad} inverted, {nt} node ties')
    blk = '\t\t<Nodes name="apex_contour_nodes">\n' + ''.join(
        f'\t\t\t<node id="{nid + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n' for i, p in enumerate(new_nodes)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    i = t.index('\t</Mesh>')
    t = t[:i] + '\t\t<Elements type="hex8" name="apex_contour">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hexes)) + '\t\t</Elements>\n' + t[i:]
    mat = re.search(r'<SolidDomain name="canal_wrap_hex_avw" mat="([^"]+)"', t).group(1)
    i = t.index('\t</MeshDomains>')
    t = t[:i] + f'\t\t<SolidDomain name="apex_contour" mat="{mat}" />\n' + t[i:]
    cb = re.search(r'<body_load name="wall_damping_canal_wrap_hex_avw"[^>]*>\s*<C>([^<]+)</C>', t)
    if cb:
        i = t.index('\t</Loads>')
        t = t[:i] + (f'\t\t<body_load name="wall_damping_apex_contour" type="mass damping" elem_set="apex_contour">\n'
                     f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n') + t[i:]
    cx = ('\t\t<constraint name="apex_contour_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
          '\t\t\t<maxaug>10</maxaug>\n' + ''.join(
              f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{c}">1</node>\n'
              + ''.join(f'\t\t\t\t<node id="{w}" bc="{c}">{-cw:.6g}</node>\n' for w, cw in ws)
              + '\t\t\t</linear_constraint>\n' for m, ws in cons for c in 'xyz') + '\t\t</constraint>\n')
    i = t.index('\t</Constraints>')
    t = t[:i] + cx + t[i:]
    msg = ("NOT IN SOURCE (the user, 2026-10-09: close the tube's apical corners; shape approved from the preview "
           "claude_diag/apex_view/apex_proposal2.html; tools/apex_contour.py): the curves' U carried on up the canal, each layer "
           "wrapped round the cervix / PVW lateral profile at its own depth, narrowing with the cervix to its tip, the U's mouth "
           'closed over the last 4 layers: ' + '; '.join(log) + f'; {len(new_nodes)} nodes, {len(hexes)} hexes in apex_contour '
           f'({mat} + its damping), {len(cons)} node ties apex_contour_ties (penalty 10, as apex_arcs)')
    return t, msg


if __name__ == '__main__':
    if sys.argv[1] == '--preview':
        preview(sys.argv[2], sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else None)
    else:
        from build_batch162 import write
        src, name = sys.argv[1], sys.argv[2]
        t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
        t2, msg = build(t)
        print(msg)
        if '--dry' not in sys.argv:
            write(name, t2, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')
