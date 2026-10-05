"""Batch 114 (2026-09-29): the PM as a structure carrying load. L113_springs_pm_step2 (the PM a clamped soft shell, the
PM / PM-middle / PM_PeB / PM_avw_bottom anchors tied to it) + the distal anterior support of batch 108 at x1e4 (the
PM_avw_bottom connectors at 1/10 of the PM curve; NOT IN SOURCE as given). One change from L113_springs_pm_step2; the
fixed-anchor counterpart is L108_springs_newline_rhoi0_avwbot1e4 (Ba -11.7 mm vs -3.3). How much of the support's lift
survives when its anchors hang on the soft PM?
  L114_springs_pm_step2_avwbot1e4
usage: py -3.10 build_batch114.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch108 import avwbot  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
name, base = 'L114_springs_pm_step2_avwbot1e4', 'L113_springs_pm_step2'
if __name__ == '__main__':
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model7(os.path.join(RUNS, base, base + '.feb'))
    avwbot(1e4)(m)
    emit(name, m)
