"""Batch 72 (2026-09-26 ~17:20): the paper's recipe on the lofts line, held past t = 1 (the user's choice (a)).

The user, viewing L63_lofts_pc10_cu30_avwparaconn_antside0_hold (+78 deg with the loads held), asked whether keeping the
6 anterior side sphincter connectors had been tested. Only to t = 1: L53_lofts_pc10_cu30_avwparaconn (them kept) gave
+41.3 deg at t = 1 (peak +45.7 at t 0.92), L59_lofts_pc10_cu30_avwparaconn_antside0 (them cut) +73.4. This run is the
hold partner of L63_..._hold: L53_lofts_pc10_cu30_avwparaconn run on to t = 1.5 with the loads held and the settle
damping (hold_to, as L61 and L63), i.e. L63_..._hold with the 6 anterior side sphincter connectors kept.
  L65_lofts_pc10_cu30_avwparaconn_hold
usage: py -3.10 build_batch72.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch68 import hold_to

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L65_lofts_pc10_cu30_avwparaconn_hold', 'L53_lofts_pc10_cu30_avwparaconn', hold_to),)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
