"""Batch 62 (2026-09-26 ~10:50): the paper's Figure 3 recipe under the paper's conditions, on both lines.

Early readings: the source Load-LA (0.014, as the paper) turns the body faster (L54_springs_pc10_cu30_LAfull +25.9 deg at
t 0.34 vs +16.2 at t 0.30 with Load-LA 1/3). One change each:
  L55_springs_pc10_cu30_m1_LAfull        L54_springs_pc10_cu30_m1 + Load-LA 0.014 (source mass and source Load-LA: the
                                         paper's conditions on the recipe)
  L55_lofts_pc10_cu30_avwparaconn_LAfull  L53_lofts_pc10_cu30_avwparaconn + Load-LA 0.014
  L55_lofts_pc10_cu30_avwparaconn_m1      L53_lofts_pc10_cu30_avwparaconn + the source chain mass
usage: py -3.10 build_batch62.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch57 import mass1
from build_batch61 import la_full

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L55_springs_pc10_cu30_m1_LAfull', 'L54_springs_pc10_cu30_m1', la_full),
          ('L55_lofts_pc10_cu30_avwparaconn_LAfull', 'L53_lofts_pc10_cu30_avwparaconn', la_full),
          ('L55_lofts_pc10_cu30_avwparaconn_m1', 'L53_lofts_pc10_cu30_avwparaconn', mass1))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
