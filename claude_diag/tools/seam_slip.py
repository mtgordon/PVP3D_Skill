"""How far the canal's lateral seam slides (2026-09-27 night): for the 116 lateral-seam AVW / cervix edge nodes of the canal
contact (build_batch90.seam_nodes), the displacement of each node relative to the PVW / PeB surface point it faced at rest
(the point the seam tie holds it to): u_node - sum_k w_k u_k, split into the part along the PVW surface's normal (opening /
closing) and the part in its plane (sliding). In a run with the tie it is ~0; in an untied run it is how far the walls
slide past each other at their edges.
usage: py -3.10 seam_slip.py RUN [RUN ...] [--t 0.76] [--seam-from RUN0] [--split]
       --split: for the 10 % most-open seam nodes, how far each wall moved along the opening direction
       --seam-from: take the seam nodes from RUN0's .feb (same mesh; needed for runs whose canal surface was cut,
       e.g. the sliding seam of batch 98, where re-finding the edge would give the wrong nodes)
"""
import os
import sys

import numpy as np

from variants7 import Model7, JOBS
from build_batch90 import seam_nodes
sys.path.insert(0, os.path.join(JOBS, 'claude_diag', 'seam_tie'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from seam_mini import project  # noqa: E402
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
args = sys.argv[1:]
tw = None
if '--t' in args:
    i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]
split = '--split' in args
args = [a for a in args if a != '--split']
src = None
if '--seam-from' in args:
    i = args.index('--seam-from'); src = args[i + 1]; args = args[:i] + args[i + 2:]
for run in args:
    m = Model7(os.path.join(RUNS, src or run, (src or run) + '.feb'))
    seam, fP, X = seam_nodes(m)
    cen = np.array([np.mean([X[n] for n in f], axis=0) for f in fP])
    pairs = []
    for v in seam:
        near = [fP[i] for i in np.argsort(np.linalg.norm(cen - X[v], axis=1))[:12]]
        d, f, w = project(X[v], X, near)
        P = np.array([X[n] for n in f])
        nrm = np.cross(P[2] - P[0], P[3] - P[1]); nrm /= np.linalg.norm(nrm)
        pairs.append((v, f, w, nrm))
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(RUNS, run, run + '.log')) or list(range(len(x.states)))
    k = conv[-1] if tw is None else min(conv, key=lambda j: abs(x.states[j][0] - tw))
    u = x.var(k, 'displacement')
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    rel = np.array([u[idx[v]] - sum(wk * u[idx[n]] for n, wk in zip(f, w)) for v, f, w, _ in pairs])
    nn = np.array([abs(r @ p[3]) for r, p in zip(rel, pairs)])
    tt = np.array([np.linalg.norm(r - (r @ p[3]) * p[3]) for r, p in zip(rel, pairs)])
    q = lambda a: f'{np.median(a):.2f} / {np.percentile(a, 90):.2f} / {a.max():.2f}'
    print(f'{run} t {x.states[k][0]:.3f}: seam slide in plane median / p90 / max {q(tt)} mm; '
          f'normal (opening) {q(nn)} mm; {len(pairs)} nodes')
    if split:
        # which wall moves: for the 10 % of seam nodes that open most, each side's displacement along the opening
        # direction (+ = toward where the AVW / cervix node went) and the total, medians
        top = np.argsort(-nn)[:max(1, len(nn) // 10)]
        ua, up, dirs = [], [], []
        for i in top:
            v, f, w, nrm = pairs[i]
            e = rel[i] / np.linalg.norm(rel[i])
            a, p = u[idx[v]], sum(wk * u[idx[n]] for n, wk in zip(f, w))
            ua.append(a @ e); up.append(p @ e); dirs.append(e)
        print(f'   the {len(top)} most-open seam nodes, displacement along the opening direction (median): AVW / cervix '
              f'node {np.median(ua):+.1f} mm, the PVW point it faced {np.median(up):+.1f} mm; opening direction '
              f'(mean) {np.round(np.mean(dirs, axis=0), 2)}; nodes {[pairs[i][0] for i in top[:6]]}')
