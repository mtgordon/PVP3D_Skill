"""Perineal body (PeB, domain _PickedSet66) displacement at a time: mean / max |u| and mean ux, uy, uz (mm), plus the
rotation of its best-fit rigid motion (deg). usage: py -3.10 peb_move.py RUN [RUN ...] [--t 1]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
args = sys.argv[1:]
tw = 1.0
if '--t' in args:
    i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]
for run in args:
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    times = np.array([s[0] for s in x.states])
    k = int(np.argmin(np.abs(times - tw)))
    dom = [d for d in x.domains if d.get('name') == '_PickedSet66']
    if not dom:
        dom = [x.domains[1]]   # fallback: the PeB is the 2nd solid domain in these models
    nodes = np.unique(np.asarray(dom[0]['conn']).ravel())
    X = np.asarray(x.X, dtype=float)[nodes]
    u = np.asarray(x.var(k, 'displacement'), dtype=float)[nodes]
    m = np.linalg.norm(u, axis=1)
    rot = float('nan')
    if X is not None:
        P, Q = X - X.mean(0), (X + u) - (X + u).mean(0)
        U, _, Vt = np.linalg.svd(P.T @ Q)
        R = (U @ Vt).T
        rot = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    print(f'{run} t {times[k]:.3f}: PeB {len(nodes)} nodes, |u| mean {m.mean():.1f} max {m.max():.1f}, '
          f'mean u ({u[:, 0].mean():.1f}, {u[:, 1].mean():.1f}, {u[:, 2].mean():.1f}), rotation {rot:.1f} deg')
