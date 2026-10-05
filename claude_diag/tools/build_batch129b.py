"""Batch 129b (2026-09-30, the user: "scale the strength by the average rather than doing it for each spring for now").
L129_tube_r30s_pm_conninner (the 26 PM_conn springs attached to the PM's inner arc, shortened from 14.0 mm to 7.4-10.3 mm)
with PM_conn_mat's elongation axis scaled by the springs' mean new length / their old length (14.0 mm), so on average a
spring makes the same force at the same fractional stretch as before (one factor for all 26, not per spring).
  L129_tube_r30s_pm_conninner_avg   L129_tube_r30s_pm_conninner + PM_conn_mat elongation x (mean L / 14.0)
usage: py -3.10 build_batch129b.py"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from febmodel import Feb  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC, OLD = 'L129_tube_r30s_pm_conninner', 'L128_tube_r30s_pm'
NAME = SRC + '_avg'
# the PM_conn family: the 26 PM_conn springs + the middle one (PM-LA-x%stiff_mat), same source law (added 2026-10-01)
SETS = ('PM_conn', 'PM-LA-x%stiff_mat')
MATS = ('PM_conn_mat', 'PM-LA-x%stiff_mat')


def lengths(run):
    p = os.path.join(RUNS, run, run + '.feb')
    f, t = Feb(p), open(p, encoding='utf-8').read()
    out = []
    for nm in SETS:
        i = t.index(f'<DiscreteSet name="{nm}">')
        blk = t[i:t.index('</DiscreteSet>', i)]
        out += [np.linalg.norm(np.array(f.nodes[int(a)]) - np.array(f.nodes[int(b)]))
                for a, b in re.findall(r'<delem>(\d+),(\d+)</delem>', blk)]
    return np.array(out)


if __name__ == '__main__':
    d = os.path.join(RUNS, NAME)
    assert not os.path.exists(d), f'{NAME} exists; not overwriting'
    L1, L0 = lengths(SRC), lengths(OLD)
    fac = L1.mean() / L0.mean()
    p = os.path.join(RUNS, SRC, SRC + '.feb')
    t = open(p, encoding='utf-8').read()
    for mat in MATS:
        i = t.index(f'name="{mat}" type="nonlinear spring">')
        j = t.index('</discrete_material>', i)
        blk = t[i:j]
        pts = re.findall(r'<pt>([^,<]+),([^<]+)</pt>', blk)
        assert len(pts) == 11
        new = re.sub(r'<pt>([^,<]+),([^<]+)</pt>', lambda mo: f'<pt>{float(mo.group(1)) * fac:.6g},{mo.group(2)}</pt>', blk)
        t = t[:i] + new + t[j:]
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t)
    msg = (f'{NAME}: {SRC} + PM_conn_mat elongation axis x {fac:.4f} (the 26 PM_conn springs\' mean length {L1.mean():.2f} mm '
           f'after the move to the inner arc / {L0.mean():.2f} mm before; range {L1.min():.1f}-{L1.max():.1f}): the same '
           f'force at the same mean fractional stretch as the source connectors; NOT IN SOURCE\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
    print(re.findall(r'<pt>[^<]+</pt>', new)[:4])
