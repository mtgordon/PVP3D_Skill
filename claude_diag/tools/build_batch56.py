"""Batch 56 (2026-09-26 ~10:12): matched candidates on both lines (the lofts line with the AVW-Para lofts as connectors).

Batch 55 on L43_springs_parcus10_clusl50: AVW-Para 50 % and PM 50 % change nothing (+41.4, +41.0 at t 0.66 vs +41.0); side
sphincters 50 % add ~3 deg (+43.8). The lofts line with the AVW-Para connectors tracks the springs line (base +3.7 / +3.4,
P-arcus 10 % +29.5 / +31.4 at t = 1). One change each, NOT IN SOURCE (connective tissue):
  L49_lofts_pc10cu50_sphside50_avwparaconn  L47_lofts_parcus10_avwparaconn_clusl50 + side sphincters 50 % (the twin of
                                            L48_springs_pc10cu50_sphside50)
  L49_springs_parcus25_clusl50              L38_springs_parcus25 + CL/USL connectors 50 % (a milder pair)
  L49_lofts_parcus25_avwparaconn_clusl50    L44_lofts_parcus25_avwparaconn + the CL/USL lofts 50 % (its twin)
usage: py -3.10 build_batch56.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale
from build_batch50 import clusl50
from build_batch52 import clusl_lofts50

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CT = 'NOT IN SOURCE (connective tissue)'
BUILDS = (('L49_lofts_pc10cu50_sphside50_avwparaconn', 'L47_lofts_parcus10_avwparaconn_clusl50',
           lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.5, CT)),
          ('L49_springs_parcus25_clusl50', 'L38_springs_parcus25', clusl50),
          ('L49_lofts_parcus25_avwparaconn_clusl50', 'L44_lofts_parcus25_avwparaconn', clusl_lofts50))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
