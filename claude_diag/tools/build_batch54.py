"""Batch 54 (2026-09-26 ~09:30): the lofts line's twins of the springs candidates, with the AVW-Para lofts as connectors.
One change each:
  L47_lofts_parcus10_avwparaconn_clusl50  L46_lofts_parcus10_avwparaconn + the CL/USL lofts at 50 % (c1, k x 0.5): the twin
                                          of L43_springs_parcus10_clusl50                NOT IN SOURCE (connective tissue)
  L47_lofts_noparcus_avwparaconn          L46_lofts_avwparaconn + every P-arcus connector cut: can the lofts line with the
                                          AVW-Para connectors take the full cut (L36_lofts_noparcus ran away at t 0.40)?
                                          DIAGNOSTIC
usage: py -3.10 build_batch54.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale
from build_batch52 import clusl_lofts50

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L47_lofts_parcus10_avwparaconn_clusl50', 'L46_lofts_parcus10_avwparaconn', clusl_lofts50),
          ('L47_lofts_noparcus_avwparaconn', 'L46_lofts_avwparaconn', lambda m: spring_scale(m, 'Parcus_conn', 0.0)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
