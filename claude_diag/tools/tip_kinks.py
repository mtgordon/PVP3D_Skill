"""Curve-tip edge kinks (2026-10-19 night; the metric used for L337 / L341 by hand, kept as a tool). For RUN at times T: each
inner node of the right curve-tip edge row 27584-27588 (27585 / 27586 / 27587): how far it sits off the line through its two
row neighbours, change from rest (mm, |offset| now - |offset| at rest), and the row's edge spread (row
length / rest length).
usage: py -3.10 tip_kinks.py RUN [RUN ...] --t 0.3,0.5"""
import os
import sys

import numpy as np

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
ROW = [27584, 27585, 27586, 27587, 27588]


def off(P, a, b, c):
    """signed-free distance of b from the line a-c, as a vector"""
    d = (c - a) / np.linalg.norm(c - a)
    v = b - a
    return v - np.dot(v, d) * d


def main(runs, ts):
    for r in runs:
        x = Xplt(os.path.join(RUNS, r, r + '.xplt'))
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        times = np.array([s[0] for s in x.states])
        X0 = x.X
        for T in ts:
            i = int(np.argmin(abs(times - T)))
            if abs(times[i] - T) > 0.01:
                print(f'{r} t {T}: not reached (last {times.max():.3f})'); continue
            P = X0 + x.var(i, 'displacement')
            out = []
            for k in range(1, 4):
                a, b, c = (idx[ROW[k - 1]], idx[ROW[k]], idx[ROW[k + 1]])
                o0 = off(X0, X0[a], X0[b], X0[c]); o1 = off(P, P[a], P[b], P[c])
                out.append(f'{ROW[k]} {np.linalg.norm(o1) - np.linalg.norm(o0):+.2f}')
            L0 = sum(np.linalg.norm(X0[idx[ROW[k + 1]]] - X0[idx[ROW[k]]]) for k in range(4))
            L1 = sum(np.linalg.norm(P[idx[ROW[k + 1]]] - P[idx[ROW[k]]]) for k in range(4))
            print(f'{r} t {times[i]:.3f}: kinks ' + ' / '.join(out) + f' mm; edge {L1 / L0:.2f}x')


if __name__ == '__main__':
    a = sys.argv[1:]
    ts = [float(v) for v in a[a.index('--t') + 1].split(',')] if '--t' in a else [0.5]
    main([v for v in a if not v.startswith('--') and (a.index(v) == 0 or a[a.index(v) - 1] != '--t')], ts)
