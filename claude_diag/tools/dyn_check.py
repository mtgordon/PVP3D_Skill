"""Quasi-static check of a (dynamic) run + LA bulge, optionally against a static run at matching times.

usage: py -3.10 dyn_check.py RUN [STATIC_RUN]
Per selected converged time: max |u| on LA nodes / on the rest, and the fastest nodal speed |du/dt| from the
two converged states around it (a dynamic run that ends quasi-static should have small speeds at t = 1).
With STATIC_RUN: max |u_dyn - u_static| at the latest converged time both runs share (LA and rest).
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from xplt import Xplt
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

RUNS = RUNS_DIR


def load(run):
    root = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot()
    names = [d.get('name') for d in root.find('MeshDomains')]
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    la = set()   # xplt domain connectivity holds 0-based node indices, not node ids
    for i, d in enumerate(x.domains):
        if i < len(names) and names[i].startswith('LA_'):
            la.update(int(n) for c in d['conn'] for n in c)
    la_i = np.array(sorted(la))
    rest = np.setdiff1d(np.arange(len(x.node_ids)), la_i)
    conv = [i for i, s in enumerate(x.states) if s[1] == 0] or list(range(len(x.states)))
    return x, la_i, rest, conv


def main():
    run = sys.argv[1]
    x, la_i, rest, conv = load(run)
    times = np.array([x.states[s][0] for s in conv])
    print(f'{run}: {len(conv)} converged states, t = {times[0]:.4f} .. {times[-1]:.4f}')
    print(f'{"t":>7} {"max|u| LA":>10} {"max|u| rest":>12} {"max speed LA":>13} {"max speed rest":>15}  [mm, mm/s]')
    for target in (0.25, 0.29, 0.33, 0.5, 0.65, 0.8, 0.9, 1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.0):
        if target > times[-1] + 1e-9:
            break
        k = int(np.argmin(np.abs(times - target)))
        k = max(k, 1)
        u1, u0 = x.var(conv[k], 'displacement'), x.var(conv[k - 1], 'displacement')
        v = np.linalg.norm(u1 - u0, axis=1) / max(times[k] - times[k - 1], 1e-12)
        un = np.linalg.norm(u1, axis=1)
        print(f'{times[k]:7.4f} {un[la_i].max():10.2f} {un[rest].max():12.2f} {v[la_i].max():13.2f} {v[rest].max():15.2f}')
    if len(sys.argv) > 2:
        xs, _, _, cs = load(sys.argv[2])
        ts = np.array([xs.states[s][0] for s in cs])
        common = min(ts[-1], times[-1])
        kd, ks = int(np.argmin(np.abs(times - common))), int(np.argmin(np.abs(ts - common)))
        du = np.linalg.norm(x.var(conv[kd], 'displacement') - xs.var(cs[ks], 'displacement'), axis=1)
        print(f'vs {sys.argv[2]} at t = {times[kd]:.4f} / {ts[ks]:.4f}: max |u_dyn - u_static| '
              f'LA {du[la_i].max():.3f} mm, rest {du[rest].max():.3f} mm')


if __name__ == '__main__':
    main()
