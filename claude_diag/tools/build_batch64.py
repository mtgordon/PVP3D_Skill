"""Batch 64 (2026-09-26 ~11:16): toward 90 deg with the paper's load on the levator.

L54_springs_pc10_cu30_LAfull (the paper's recipe, P-arcus 10 % + CL/USL 30 %, with the source Load-LA 0.014) turns the
body +51.5 deg at t 0.68 (peak +53.1), the most of any run that keeps going; with Load-LA 1/3 the recipe reaches +44.4.
One change each (NOT IN SOURCE, connective tissue):
  L57_springs_pc0_cu30_LAfull           L56_springs_pc0_cu30 (P-arcus cut, CL/USL 30 %) + Load-LA 0.014
  L57_springs_pc10_cu30_LAfull_antside0  L54_springs_pc10_cu30_LAfull + the 3 anterior pairs of side sphincters cut
usage: py -3.10 build_batch64.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch59 import antside0
from build_batch61 import la_full

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L57_springs_pc0_cu30_LAfull', 'L56_springs_pc0_cu30', la_full),
          ('L57_springs_pc10_cu30_LAfull_antside0', 'L54_springs_pc10_cu30_LAfull', antside0))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
