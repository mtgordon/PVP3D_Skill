"""Perineal-body rotation comparison table (batch 44 on): convergence plus the body's best-fit rotation about -x.

usage: py -3.10 peb_table.py RUN [RUN ...] [--every 4]
Columns: change (the .feb.changes.txt lines after 'base:'), t reached, converged steps, failed attempts, status, wall;
rot_-x (deg, + = the user's sense: cranial end forward and down) at t 0.5 / 0.66 / 1.0 (nearest converged state, if
reached), the peak over the run with its t, and at the last converged state (with, in brackets, the turn of the body's
long principal axis in the y-z plane, the same sense: a check on the rigid fit when the body bends); the centroid
displacement (y, z) there.
--every N: sample every N-th converged state for the peak (default 4; the last state is always included).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from xplt import Xplt  # noqa: E402
from batch18_table import log_info  # noqa: E402
from peb_rotation import kabsch, region_nodes  # noqa: E402


def change_text(run):
    p = os.path.join(RUNS_DIR, run, run + '.feb.changes.txt')
    if not os.path.exists(p):
        return ''
    lines = [l.strip() for l in open(p, encoding='latin-1') if l.strip() and not l.startswith('base:')]
    return '; '.join(lines)


def axis_angle(Q):
    """Angle (deg) above the +y axis of the body's long principal axis in the y-z plane (oriented toward +y)."""
    c = Q - Q.mean(axis=0)
    w, v = np.linalg.eigh(np.cov(c.T))
    a = v[:, int(np.argmax(w))]
    if a[1] < 0:
        a = -a
    return np.degrees(np.arctan2(a[2], a[1]))


def rot_series(run, every=4):
    x = Xplt(os.path.join(RUNS_DIR, run, run + '.xplt'))
    nodes = region_nodes(run, x, ('_PickedSet66',))
    P = x.X[nodes]
    conv = [k for k, st in enumerate(x.states) if st[1] == 0]
    if not conv:
        return None
    ks = sorted(set(conv[::every]) | {conv[-1]})
    tt, rr, cu = [], [], []
    a0 = axis_angle(P)
    for k in ks:
        u = x.var(k, 'displacement')[nodes]
        _, rv, _ = kabsch(P, P + u)
        tt.append(x.states[k][0]); rr.append(-np.degrees(rv[0])); cu.append(u.mean(axis=0))
    ax_end = a0 - axis_angle(P + u)
    return np.array(tt), np.array(rr), np.array(cu), ax_end


def main():
    args = sys.argv[1:]
    every = 4
    if '--every' in args:
        i = args.index('--every')
        every = int(args[i + 1])
        args = args[:i] + args[i + 2:]
    print('| run | change | t reached | steps | failed | status | wall | rot_-x t 0.5 / 0.66 / 1.0 | peak (t) | '
          'at end | centroid u y, z at end |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for run in args:
        try:
            t, n, fl, nj, st, wall = log_info(run)
        except FileNotFoundError:
            print(f'| {run} | {change_text(run)} | (no log) | | | | | | | | |')
            continue
        s = rot_series(run, every)
        if s is None:
            print(f'| {run} | {change_text(run)} | {t:.4f} | {n} | {fl} | {st} | {wall} | | | | |')
            continue
        tt, rr, cu, ax_end = s

        def at(tq):
            if tq > tt[-1] + 1e-3:
                return '-'
            return f'{rr[int(np.argmin(np.abs(tt - tq)))]:+.1f}'
        ip = int(np.argmax(rr))
        print(f'| {run} | {change_text(run)} | {t:.4f} | {n} | {fl} | {st} | {wall} | {at(0.5)} / {at(0.66)} / '
              f'{at(1.0)} | {rr[ip]:+.1f} ({tt[ip]:.2f}) | {rr[-1]:+.1f} (axis {ax_end:+.1f}) | {cu[-1][1]:+.1f}, {cu[-1][2]:+.1f} |')


if __name__ == '__main__':
    main()
