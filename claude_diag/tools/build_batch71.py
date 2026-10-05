"""Batch 71 (2026-09-26 ~13:40): where the best set settles on the source chain mass (the most paper-like runs).

L63_lofts_pc10_cu30_avwparaconn_antside0_m1 turns +80.0 deg at t 0.71 (8 failed); L60_springs_pc10_cu30_antside0_m1 peaks
+73.5 at t 0.85; the hold run L61_springs_pc10_cu30_antside0_hold shows the t = 1 state is a low point of a swing
(+67.2 at t = 1, +73.3 at t 1.10). One change each (the run length, as build_batch68.hold_to):
  L64_lofts_pc10_cu30_avwparaconn_antside0_m1_hold  L63_lofts_pc10_cu30_avwparaconn_antside0_m1 run on to t = 1.5
  L64_springs_pc10_cu30_antside0_m1_hold            L60_springs_pc10_cu30_antside0_m1 run on to t = 1.5
usage: py -3.10 build_batch71.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch68 import hold_to

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L64_lofts_pc10_cu30_avwparaconn_antside0_m1_hold', 'L63_lofts_pc10_cu30_avwparaconn_antside0_m1'),
          ('L64_springs_pc10_cu30_antside0_m1_hold', 'L60_springs_pc10_cu30_antside0_m1'))
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        hold_to(m)
        emit(name, m)
