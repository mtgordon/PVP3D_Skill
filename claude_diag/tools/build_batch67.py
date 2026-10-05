"""Batch 67 (2026-09-26 ~12:05): the anterior side sphincter cut is the smooth way past 60 deg.

L57_springs_pc10_cu30_LAfull_antside0 (the paper's recipe + Load-LA 0.014 + the 3 anterior pairs of side sphincters cut)
peaks +74.7 deg at t 0.58 with 4 failed attempts; L56_springs_pc10_cu30_antside0 (Load-LA 1/3) +71.0 at t 0.87 with 7.
The P-arcus-cut and full-Load-LA runs without it snap past ~55-66 deg (a 3 s ramp does not help). One change each:
  L60_lofts_pc10_cu30_avwparaconn_antside0_LAfull  L59_lofts_pc10_cu30_avwparaconn_antside0 + Load-LA 0.014 (the lofts
                                                   twin of L57_springs_pc10_cu30_LAfull_antside0)
  L60_springs_pc10_cu30_antside0_m1               L56_springs_pc10_cu30_antside0 + the source chain mass
  L60_springs_pc10_cu30_LAfull_antside0_m1        L57_springs_pc10_cu30_LAfull_antside0 + the source chain mass (the
                                                   paper's conditions on the best set)
usage: py -3.10 build_batch67.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch57 import mass1
from build_batch61 import la_full

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L60_lofts_pc10_cu30_avwparaconn_antside0_LAfull', 'L59_lofts_pc10_cu30_avwparaconn_antside0', la_full),
          ('L60_springs_pc10_cu30_antside0_m1', 'L56_springs_pc10_cu30_antside0', mass1),
          ('L60_springs_pc10_cu30_LAfull_antside0_m1', 'L57_springs_pc10_cu30_LAfull_antside0', mass1))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
