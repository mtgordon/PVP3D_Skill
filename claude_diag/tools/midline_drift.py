"""Lateral (x) drift of the vaginal-wall midline (Abaqus XSYMM set BC-VW-mid) over converged states.
usage: py -3.10 midline_drift.py run1 [run2 ...]"""
import os, re, sys
import numpy as np
from xplt import Xplt
from febmodel import Feb
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

R = RUNS_DIR
for run in sys.argv[1:]:
    d = os.path.join(R, run)
    f = Feb(os.path.join(d, run + '.feb'))
    x = Xplt(os.path.join(d, run + '.xplt'))
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    mid = [idx[n] for n in f.nodesets['BC-VW-mid'] if n in idx]
    log = open(os.path.join(d, run + '.log')).read()
    ct = [float(v) for v in re.findall(r'^------- converged at time : ([0-9.e-]+)', log, re.M)]
    times = np.array([s[0] for s in x.states])
    print(f'{run}: {len(ct)} converged steps')
    for t in sorted(set(ct[::max(1, len(ct) // 8)] + ct[-2:])):
        js = np.where(np.abs(times - t) < 1e-6 * max(t, 1e-3))[0]
        if not len(js):
            continue
        u = x.var(int(js.max()), 'displacement')
        um = u[mid]
        print(f'   t={t:.4f}  midline |ux| max {np.abs(um[:, 0]).max():.4f}  mean ux {um[:, 0].mean():+.4f}   '
              f'|u| max {np.linalg.norm(um, axis=1).max():.3f}   whole-model |ux| max {np.abs(u[:, 0]).max():.3f}')
