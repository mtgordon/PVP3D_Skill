"""Batch 94 (2026-09-27 ~23:40): the canal seam tie (NOT IN SOURCE; the user's to adopt or not) carried on.

L86_newline_vwyeoh_rhoi0_seamtie (the new springs line + the Yeoh walls + rhoi 0 + the lateral seam tie) went from t 0.785
to 0.990 in 20 min with 3 failed attempts, while every contact / solver setting crawled or failed past t 0.78-0.89.
Carried on, one change each from an existing run:
  L89_springs_newline_rhoi0_seamtie        L87_springs_newline_rhoi0 (the new springs line + rhoi 0, Ogden walls) + the tie:
                                           how much the tie itself changes the line's solution
  L89_newline_vwyeoh_seamtie               L84_springs_newline_vwyeoh (the faithful walls, rhoi 0.5) + the tie: is rhoi 0
                                           still needed once the seam is tied?
  L89_lofts_newline_vwyeoh_rhoi0_seamtie   L87_lofts_newline_vwyeoh_rhoi0 (the new lofts line + the faithful walls + rhoi 0,
                                           crawled at t 0.52) + the tie
usage: py -3.10 build_batch94.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch90 import seam_tie

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

BUILDS = (('L89_springs_newline_rhoi0_seamtie', 'L87_springs_newline_rhoi0', seam_tie),
          ('L89_newline_vwyeoh_seamtie', 'L84_springs_newline_vwyeoh', seam_tie),
          ('L89_lofts_newline_vwyeoh_rhoi0_seamtie', 'L87_lofts_newline_vwyeoh_rhoi0', seam_tie))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        print(name, ':', m.log[-1][:150])
        emit(name, m)
