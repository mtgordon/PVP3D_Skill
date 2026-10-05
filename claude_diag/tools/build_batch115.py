"""Batch 115 (2026-09-29): batch 110's edge holds rebuilt after the bug fix in build_batch109._hold_edges (a run of
free edge nodes that goes round more than 25 % of the loop is no longer filled: in batch 110 / 112 it fixed USL-L_fan's
whole boundary). Holds as batch 110 (augmented linear constraints, penalty 10, maxaug 10). NOT IN SOURCE.
  L115_fibre_all_m05_edges   L109_fibre_all_m05 + every loft's edges held (the user's choice, fibre lofts with held edges)
  L115_fibre_usl_m05_edges   L109_fibre_usl_m05 + the USL lofts' edges held
usage: py -3.10 build_batch115.py [NAME ...]"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch109 import Model109, fibre, ALL, BASE, RUNS, emit  # noqa: E402

WANT = set(sys.argv[1:])
BUILDS = [('L115_fibre_all_m05_edges', BASE, fibre(ALL, 0.05, edges=True)),
          ('L115_fibre_usl_m05_edges', BASE, fibre(('usl',), 0.05, edges=True))]
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        if os.path.exists(os.path.join(RUNS, name)):
            print(f'{name} exists; not overwriting')
            continue
        mdl = Model109(os.path.join(RUNS, base, base + '.feb'))
        fn(mdl)
        emit(name, mdl)
