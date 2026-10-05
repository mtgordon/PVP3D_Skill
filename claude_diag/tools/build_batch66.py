"""Batch 66 (2026-09-26 ~11:50): the anterior side sphincter cut on the paper's recipe, on both lines.

L56_springs_pc10_cu30_antside0 (the paper's recipe, P-arcus 10 % + CL/USL 30 %, plus the 3 anterior pairs of side
sphincter connectors cut) turns the body +58.7 deg at t 0.5 and +66.7 at t 0.66 (peak +67.5) with 3 failed attempts: the
best smooth run. One change each, NOT IN SOURCE (connective tissue):
  L59_lofts_pc10_cu30_avwparaconn_antside0  L53_lofts_pc10_cu30_avwparaconn + the anterior side sphincters cut (its twin)
  L59_springs_pc0_cu30_antside0             L56_springs_pc0_cu30 (P-arcus cut) + the same
usage: py -3.10 build_batch66.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch59 import antside0

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L59_lofts_pc10_cu30_avwparaconn_antside0', 'L53_lofts_pc10_cu30_avwparaconn'),
          ('L59_springs_pc0_cu30_antside0', 'L56_springs_pc0_cu30'))
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        antside0(m)
        emit(name, m)
