"""Batch 130 (2026-09-30 night, the user: "If they work, do some impaired runs"): the paper's impaired cases P1 / P2 (Luo
et al. 2015; build_batch107.py: levator 80 / 40 % of healthy, apical CL / USL x 0.7 / 0.4, Parcus x 0.15) on the new
base L128_tube_r30s_pm (the tube + the PM as a structure) and on its two PM_conn variants (the PM_conn springs attached to
the PM's inner arc, with PM_conn_mat as it was or scaled by the mean length change).
  L130_tube_pm_paperP1 / P2                L128_tube_r30s_pm + P1 / P2
  L130_tube_pm_conninner_paperP1 / P2      L129_tube_r30s_pm_conninner + P1 / P2
  L130_tube_pm_conninner_avg_paperP1 / P2  L129_tube_r30s_pm_conninner_avg + P1 / P2
  L130_tube_pm_wallcontact_paperP1 / P2    L131_tube_pm_wallcontact (+ wall-LA and wall-PM contact) + P1 / P2
  L130_thin275_wrapend_paperP1 / P2        L143_thin275_narrow3_wrapend_noflare + P1 / P2 (2026-10-01)
usage: py -3.10 build_batch130.py NAME [NAME ...]"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch107 import P1, P2  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BUILDS = {}
for tag, base in (('tube_pm', 'L128_tube_r30s_pm'), ('tube_pm_conninner', 'L129_tube_r30s_pm_conninner'),
                  ('tube_pm_conninner_avg', 'L129_tube_r30s_pm_conninner_avg'),
                  ('tube_pm_wallcontact', 'L131_tube_pm_wallcontact'),
                  ('thin275_wrapend', 'L143_thin275_narrow3_wrapend_noflare')):
    BUILDS[f'L130_{tag}_paperP1'] = (base, P1)
    BUILDS[f'L130_{tag}_paperP2'] = (base, P2)

if __name__ == '__main__':
    for name in sys.argv[1:]:
        base, fn = BUILDS[name]
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
