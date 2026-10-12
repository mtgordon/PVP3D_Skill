"""Preview (no model written), 2026-10-09 after the user's "there's still a hole at the end of the U. Also, it doesn't contour
to the cervix quite as nicely as I'd like": the curves' U carried on up the canal, but every layer wraps the walls' lateral
profile AT ITS OWN DEPTH. Sections of the cervix + PVW outer surfaces give, per depth, the lateral face's position (x_w) and AP
extent (y_lo..y_hi, the 1 mm strip next to the most lateral point). Each cap node keeps its relative AP position s across that
extent and its lateral offset h from the face; the bulge shrinks with the extent (not below x0.45, so the shell keeps ~1 mm) as
the cervix rounds off to its tip, so the U narrows and closes onto the tip; over the last CLOSE layers the U's inner surface comes down onto the walls' face, so
the end is solid (no opening in the U's mouth). Output: JSON for apex_view (part 'proposal').
usage: py -3.10 apex_contour_preview.py RUN OUT.json [ZEND]"""
import json
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apex_close import HF  # noqa: E402
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


def main(run, out, zend=None):
    mesh = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot().find('Mesh')
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
    nodes, quads = [], []
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
            P0 = P0[::-1]                                   # row 0 = the U's inner surface (smaller arc)
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
        base = len(nodes)
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
                    nodes.append(p)
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
        print(f'{"right" if sgn > 0 else "left"}: cap depth {P0[:, :, 2].min():.1f}..{P0[:, :, 2].max():.1f} -> end {ze:.2f}, '
              f'AP extent {span0.mean():.1f} -> {f(ze)[2] - f(ze)[1]:.1f} mm, {n} layers')
    json.dump({'nodes': [[round(float(c), 3) for c in p] for p in nodes], 'parts': {'proposal': quads}, 'edges': {}},
              open(out, 'w'))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else None)
