"""Curve-element swelling (2026-10-10: the tube main's MAX curves balloon at their distal PVW ends, "exploded at the end").
For each run: the volume ratio V / V0 of every canal_wrap_hex_* element (hex8 split into 5 tets, from the plotted
displacements) at the last converged state (or the last one at or before --t), the max ratio and its element, how many
elements exceed 2x, the smallest ratios (a negative one = the hex is folded inside out somewhere), and the max ratio per 0.1 of time (when the swelling starts).
usage: py -3.10 curve_swell.py RUN [RUN ...] [--t 0.9] [--prefix canal_wrap_hex]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
TETS = [(0, 1, 2, 5), (0, 2, 3, 7), (0, 5, 7, 4), (2, 7, 5, 6), (0, 2, 7, 5)]  # 5-tet split of a hex8


def vols(X, conn):
    v = np.zeros(len(conn))
    for a, b, c, d in TETS:
        p = X[conn[:, a]]
        v += np.einsum('ij,ij->i', X[conn[:, b]] - p, np.cross(X[conn[:, c]] - p, X[conn[:, d]] - p)) / 6.0
    return v


def report(run, tw, prefix):
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    doms = [d for d in x.domains if str(d.get('name', '')).startswith(prefix) and np.asarray(d['conn']).shape[1] == 8]
    if not doms:
        return f'| `{run}` | no {prefix}* hex domains in the plot |'
    conn = np.vstack([np.asarray(d['conn']) for d in doms])
    eids = np.concatenate([np.asarray(d['eids']) for d in doms])
    X0 = x.X
    V0 = vols(X0, conn)
    idx = converged_state_indices(x, os.path.join(RUNS, run, run + '.log'))
    if tw is not None:
        idx = [i for i in idx if x.states[i][0] <= tw + 1e-9]
    if not idx:
        return f'| `{run}` | no converged state |'
    hist = {}
    for i in idx[::max(1, len(idx) // 40)] + [idx[-1]]:
        t = x.states[i][0]
        r = vols(X0 + x.var(i, 'displacement'), conn) / V0
        b = round(np.floor(t * 10 + 1e-9) / 10, 1)
        hist[b] = max(hist.get(b, 0), r.max())
        last = (t, r)
    t, r = last
    k = np.argsort(-r)[:6]
    worst = ', '.join(f'{eids[j]} {r[j]:.1f}' for j in k)
    folded = ', '.join(f'{eids[j]} {r[j]:.2f}' for j in np.argsort(r)[:3])
    trend = ' '.join(f'{b:.1f}:{v:.1f}' for b, v in sorted(hist.items()))
    return (f'| `{run}` | {t:.3f} | {r.max():.2f} | {(r > 2).sum()} | {(r > 1.5).sum()} | {folded} | {worst} | {trend} |')


def main(argv):
    tw, prefix, runs = None, 'canal_wrap_hex', []
    while argv:
        a = argv.pop(0)
        if a == '--t':
            tw = float(argv.pop(0))
        elif a == '--prefix':
            prefix = argv.pop(0)
        else:
            runs.append(a)
    print('| run | t | max V/V0 | > 2x | > 1.5x | min V/V0 (5-tet split; < 0 = folded) | worst (element ratio) | max V/V0 by time |')
    print('|---|---|---|---|---|---|---|---|')
    for run in runs:
        print(report(run, tw, prefix))


if __name__ == '__main__':
    main(sys.argv[1:])
