"""Batch 85 (2026-09-27 ~18:30): which PVW_LA setting for the lines.
  L81_lofts_la3_nopinch_pvwpen05       L74_lofts_la3_nopinch (the lofts line + option (b)) + PVW_LA penalty 5 -> 0.5 (no
                                       seg_up): does the lofts line keep its speed? (with seg_up 2: 36 failed, 10985
                                       iterations vs 3068 for L26_lofts_la3_pm)
  L81_springs_la3_nopinch_pvwsegup2_pen05  L77_springs_la3_nopinch_pvwsegup2 (the springs line) + PVW_LA penalty 5 -> 0.5:
                                       the Ogden-wall control of L79_springs_la3_vwyeoh_pvwpen05 (the walls the only difference)
usage: py -3.10 build_batch85.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def pvwpen05(m):
    assert m.set_contact('PVW_LA', penalty=0.5) == 1


BUILDS = (('L81_lofts_la3_nopinch_pvwpen05', 'L74_lofts_la3_nopinch', pvwpen05),
          ('L81_springs_la3_nopinch_pvwsegup2_pen05', 'L77_springs_la3_nopinch_pvwsegup2', pvwpen05))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
