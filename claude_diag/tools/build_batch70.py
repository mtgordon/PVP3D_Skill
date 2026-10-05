"""Batch 70 (2026-09-26 ~13:10): confirm the best lofts model.

L59_lofts_pc10_cu30_avwparaconn_antside0 (the lofts line with the AVW-Para connectors, P-arcus 10 %, the CL/USL lofts
30 %, the 3 anterior pairs of side sphincters cut) reached t = 1 with +73.4 deg (peak +74.8), 21 failed attempts: the
largest rotation at full load. One change each:
  L63_lofts_pc10_cu30_avwparaconn_antside0_hold  run on to t = 1.5 with the loads held and the settle damping (as L61)
  L63_lofts_pc10_cu30_avwparaconn_antside0_m1    the source arcus chain mass
usage: py -3.10 build_batch70.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch57 import mass1
from build_batch68 import hold_to

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BASE = 'L59_lofts_pc10_cu30_avwparaconn_antside0'
BUILDS = (('L63_lofts_pc10_cu30_avwparaconn_antside0_hold', hold_to),
          ('L63_lofts_pc10_cu30_avwparaconn_antside0_m1', mass1))
if __name__ == '__main__':
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, BASE, BASE + '.feb'))
        fn(m)
        emit(name, m)
