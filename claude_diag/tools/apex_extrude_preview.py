"""Preview (no model written): the curves' own U cross-section carried on unchanged up the canal to the cervix's end
(2026-10-09, after the user's "make the radius stay constant" and "keep going even though the pvw is not there on the other
side"). At each corner the curves' apical end cap (5 x 7 nodes: 5 through the thickness, 7 round the U) is swept along the
curves' own centre-line direction; every node goes on until it reaches the depth of the cervix's tip, so the swept part ends
flat at the cervix's end. Output: JSON for apex_view (part 'proposal').
usage: py -3.10 apex_extrude_preview.py RUN OUT.json [--step 0.7]"""
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


def main(run, out, step=0.7):
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
        assert Gr.shape == (5, 7), Gr.shape
        P0 = np.array([[L[v] for v in r] for r in Gr])
        d0 = P0[:, :, 2].max()
        # the curves' centre-line direction near their end: centroid of the last 2 mm vs the 2 mm before 5 mm back
        cv = np.array([L[v] for d in CURVES for el in E[d] for v in el])
        cv = cv[sgn * cv[:, 0] > 1.5]
        a = cv[(cv[:, 2] > d0 - 2.5)].mean(0); b = cv[(cv[:, 2] > d0 - 7.5) & (cv[:, 2] < d0 - 5)].mean(0)
        dirn = (a - b) / np.linalg.norm(a - b)
        cx = np.array([L[v] for v in cxn]); cx = cx[sgn * cx[:, 0] > 0]
        dtip = cx[:, 2].max()
        tlen = (dtip - P0[:, :, 2]) / dirn[2]                 # per node: the sweep length to reach the cervix tip's depth
        n = max(2, int(np.ceil(tlen.max() / step)))
        base = len(nodes)
        for k in range(n + 1):
            for i in range(5):
                for j in range(7):
                    nodes.append(P0[i, j] + dirn * tlen[i, j] * k / n)
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
        print(f'{"right" if sgn > 0 else "left"}: cap depth {P0[:, :, 2].min():.1f}..{d0:.1f}, direction '
              f'({dirn[0]:+.2f}, {dirn[1]:+.2f}, {dirn[2]:+.2f}), to the cervix tip depth {dtip:.1f}: sweep '
              f'{tlen.min():.1f}..{tlen.max():.1f} mm, {n} layers')
    json.dump({'nodes': [[round(float(c), 3) for c in p] for p in nodes], 'parts': {'proposal': quads}, 'edges': {}},
              open(out, 'w'))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
