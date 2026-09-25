"""Signed distance of PVW load-surface nodes to the deformed LA mid-surface at converged states."""
import os, re, sys
import numpy as np
from xplt import Xplt
from febmodel import Feb, pt_tri_dist
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

run = sys.argv[1]
d = os.path.join(RUNS_DIR, run)
f = Feb(os.path.join(d, run + '.feb'))
x = Xplt(os.path.join(d, run + '.xplt'))
idx = {int(n): i for i, n in enumerate(x.node_ids)}
tris = []
for name in ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri'):
    for c in f.elem_blocks[name][1].values():
        tris += [c[:3]] if len(c) == 3 else [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]
pv = [n for n in f.surface_nodes('Load-PVW_surf') if n in idx]
log = open(os.path.join(d, run + '.log')).read()
ct = [float(v) for v in re.findall(r'^------- converged at time : ([0-9.e-]+)', log, re.M)]
times = np.array([s[0] for s in x.states])
for t in ct[::max(1, len(ct) // 6)] + ct[-1:]:
    js = np.where(np.abs(times - t) < 1e-6 * max(t, 1e-3))[0]
    if not len(js):
        continue
    u = x.var(int(js.max()), 'displacement')
    X = x.X + u
    A = np.array([X[idx[a]] for a, b, c in tris]); B = np.array([X[idx[b]] for a, b, c in tris]); C = np.array([X[idx[c]] for a, b, c in tris])
    cen = (A + B + C) / 3
    dmin = []
    for n in pv:
        p = X[idx[n]]
        near = np.argsort(np.linalg.norm(cen - p, axis=1))[:30]
        dmin.append(min(pt_tri_dist(p, A[i], B[i], C[i]) for i in near))
    dmin = np.array(dmin)
    print(f't={t:.4f}: PVW-back to LA mid-surface min {dmin.min():.2f} mm, nodes < 2 mm (inside LA thickness): {int((dmin < 2).sum())}')
