"""Batch 128 (2026-09-30, the user: "make a new base model that includes that as well"): the new base L127_tube_r30s (the
tube, faithful walls, seg_up 2) + the perineal membrane as a structure, as batch 113's step 2 (build_batch113.py):
PM_Plane a deformable 2 mm shell (PM_Yeoh, the PM_mid PARAVAG_H_highdensity fit), its outer arc clamped (as OPAL325_PM_mid),
the rigid-body constraint removed, and the PM_PeB, PM_avw_bottom, PM and PM-LA anchors (47) tied to it (penalty 100,
maxaug 0), their fixed BCs removed. NOT IN SOURCE (the source's PM_Plane is a display body; the anchors are fixed points).
On the older springs line (L87_springs_newline_rhoi0) the same change left Ba / Bp / C unchanged (L113_springs_pm_step2).
  L128_tube_r30s_pm   L127_tube_r30s + the PM as a structure with its anchors tied to it
usage: py -3.10 build_batch128.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS, emit  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L127_tube_r30s'

if __name__ == '__main__':
    name = 'L128_tube_r30s_pm'
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model113(os.path.join(RUNS, BASE, BASE + '.feb'))
    step2(INNER + OUTER_BOTTOM)(m)
    for line in m.log[-3:]:
        print(line[:400])
    emit(name, m)
