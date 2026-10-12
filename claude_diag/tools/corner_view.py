"""Render the tube's apical corners (boundary faces of cervix, PVW, curves incl. *_apex) and ray-test the apex lumen for
holes. usage: py -3.10 corner_view.py RUN OUT.png"""
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apex_close import _ray_tri, HF  # noqa: E402
from build_batch202 import plane  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
run, out = sys.argv[1], sys.argv[2]
mesh = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot().find('Mesh')
X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')}
c0, nrm = plane()
nn = nrm if np.dot(np.mean([X[v] for el in E['_PickedSet346'] for v in el], 0) - c0, nrm) < 0 else -nrm
e1 = np.array([1, 0, 0.]); e1 -= e1 @ nn * nn; e1 /= np.linalg.norm(e1); e2 = np.cross(nn, e1)
L = lambda p: np.array([(p - c0) @ e1, (p - c0) @ e2, (p - c0) @ nn])
col = {'_PickedSet64': 'red', '_PickedSet346': 'green', '_PickedSet347': 'blue', '_PickedSet66': 'gray'}
for d in E:
    if d.startswith('canal_wrap_hex_avw'):
        col[d] = 'gold' if d.endswith('_apex') else 'orange'
    if d.startswith('canal_wrap_hex_pvw'):
        col[d] = 'cyan' if d.endswith('_apex') else 'purple'
cnt, own = Counter(), {}
for d in col:
    for el in E.get(d, []):
        for fc in HF:
            k = frozenset(el[i] for i in fc); cnt[k] += 1; own[k] = (d, [el[i] for i in fc])
bnd = [own[k] for k, c in cnt.items() if c == 1]
fig = plt.figure(figsize=(16, 8))
for k_, (side, az) in enumerate(((1, -60), (-1, -120))):
    ax = fig.add_subplot(1, 2, k_ + 1, projection='3d')
    polys, cols = [], []
    for d, fc in bnd:
        Q = np.array([L(X[n]) for n in fc]); m = Q.mean(0)
        if m[2] < -78 and side * m[0] > 12:
            polys.append(Q); cols.append(col[d])
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cols, edgecolors='k', linewidths=0.2, alpha=0.85))
    A = np.vstack(polys); ax.set_xlim(A[:, 0].min(), A[:, 0].max()); ax.set_ylim(A[:, 1].min(), A[:, 1].max())
    ax.set_zlim(A[:, 2].min(), A[:, 2].max()); ax.view_init(elev=-25, azim=az)
    ax.set_title(('right' if side > 0 else 'left') + ' apical corner (gold / cyan = new apex hexes)')
plt.tight_layout(); plt.savefig(out, dpi=80)
tris = []
for d, fc in bnd:
    q = [X[n] for n in fc]; tris += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
T = np.array(tris)
rng = np.random.default_rng(0)
Dd = rng.normal(size=(800, 3)); Dd /= np.linalg.norm(Dd, axis=1)[:, None]; Dd = Dd[(Dd @ nn) < -0.2]
for hgt, ey, lat in ((-80, 4.5, 0), (-80, 4.5, 16), (-85, 3.0, 16), (-85, 3.0, -16), (-70, 5.5, 16)):
    o = c0 + hgt * nn + ey * e2 + lat * e1
    print(f'lumen point height {hgt}, lateral {lat}: {sum(_ray_tri(o, d, T) is None for d in Dd)} of {len(Dd)} apical rays escape')
