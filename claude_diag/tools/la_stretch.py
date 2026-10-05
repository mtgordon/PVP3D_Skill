"""How much the LA stretches (balloons), not how far it moves: per LA element the area stretch A/A0 (median / p90 / max),
the total LA area ratio, and the element-edge stretch, at the converged states nearest the requested times.
The user's concern (2026-09-27) is this stretching; at the source loads L67_pen5_ctrl reaches a total LA area x2.16.
usage: py -3.10 la_stretch.py RUN [RUN ...] [--t 0.5,1.0]   (default t = the last converged state)"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
args = sys.argv[1:]
ts = None
if '--t' in args:
    i = args.index('--t'); ts = [float(v) for v in args[i + 1].split(',')]; args = args[:i] + args[i + 2:]


def area(P):
    if len(P) == 4:
        return 0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[1]))
    return 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))


for run in args:
    root = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot()
    els = [[int(v) for v in e.text.split(',')] for b in root.find('Mesh').findall('Elements')
           if b.get('name', '').startswith('LA_') for e in b]
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    T = np.array([s[0] for s in x.states])
    ids = [[idx[n] for n in c] for c in els]
    A0 = np.array([area(x.X[c]) for c in ids])
    L0 = np.array([np.linalg.norm(x.X[c[i]] - x.X[c[(i + 1) % len(c)]]) for c in ids for i in range(len(c))])
    for k in ([int(np.argmin(np.abs(T - t))) for t in ts] if ts else [len(T) - 1]):
        X = x.X + x.var(k, 'displacement')
        r = np.array([area(X[c]) for c in ids]) / A0
        st = np.array([np.linalg.norm(X[c[i]] - X[c[(i + 1) % len(c)]]) for c in ids for i in range(len(c))]) / L0
        print(f'{run} t={T[k]:.3f}: LA area x{(r * A0).sum() / A0.sum():.2f}; element area stretch median {np.median(r):.2f} / '
              f'p90 {np.percentile(r, 90):.2f} / max {r.max():.2f}; edge stretch median {np.median(st):.2f} / p90 '
              f'{np.percentile(st, 90):.2f}')
