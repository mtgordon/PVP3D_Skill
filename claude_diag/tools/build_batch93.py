"""Batch 93 (2026-09-27 ~23:05): canal seg_up 5, the best setting for the faithful walls so far, checked for adoption.

At 23:00 L86_newline_vwyeoh_rhoi0_se1segup5 (the faithful walls + rhoi 0 + canal seg_up 5) led at t 0.859 and still moved;
dtol 0.01 (0.829) and PVW_LA penalty 0.25 (0.786) crawled. rhoi 0 alone keeps the new springs line's solution at t = 1
(L87_springs_newline_rhoi0: Ba / Bp / C -3.3 / -19.2 / -38.3 vs -3.2 / -19.2 / -38.3, LA x1.09) and is 28 % faster. One
change each:
  L88_springs_newline_rhoi0_se1segup5        L87_springs_newline_rhoi0 + SlidingElastic1 seg_up 0 -> 5: does the line (Ogden
                                             walls) keep its solution and speed with it?
  L88_lofts_newline_vwyeoh_rhoi0_se1segup5   L87_lofts_newline_vwyeoh_rhoi0 (the lofts line + the faithful walls + rhoi 0,
                                             t 0.52 at 23:00, slow) + SlidingElastic1 seg_up 0 -> 5
usage: py -3.10 build_batch93.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch79 import se1

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

BUILDS = (('L88_springs_newline_rhoi0_se1segup5', 'L87_springs_newline_rhoi0', se1(seg_up=5)),
          ('L88_lofts_newline_vwyeoh_rhoi0_se1segup5', 'L87_lofts_newline_vwyeoh_rhoi0', se1(seg_up=5)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
