"""Batch 50 (2026-09-26 ~08:55): getting the lofts line through its crawl, and the apex on the springs line.

The lofts line crawls at t ~0.43-0.52 whenever the proximal P-arcus is weakened (parcus50 got through in 63 min; parcus25,
parcusprox0, the USL-L-merged parcus50 did not). L38_lofts_parcus25 at t 0.434: no distortion (J 0.99-1.06, |u| <= 9.6 mm);
the fastest nodes in the AVW and the AVW-Para lofts (~40 mm/s); the failed steps fail on the displacement norm with the
energy converged, 75 zero line steps: a floppy mode. AVW-Para-L is the other loft still tied by tied-node-on-facet.
One change each:
  L43_lofts_parcus25_avwparatie   L38_lofts_parcus25 + the AVW-Para-L loft on shared nodes (merge_tnof; equivalent if every
                                  tied tissue node sits on a loft node, which merge_tnof checks); a loft convention
  L43_lofts_parcus25_dtol01       L38_lofts_parcus25 + solver dtol 0.001 -> 0.01 (solver only; etol still 0.01)
  L43_springs_parcus10_clusl50    L40_springs_parcus10 + every CL and USL connector set at 50 %   NOT IN SOURCE (connective
                                  tissue): with the PVW's lateral hold weak, does a lower apex turn the body further?
usage: py -3.10 build_batch50.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale
from tnof_merge import merge_tnof

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CT = 'NOT IN SOURCE (connective tissue)'
AVW_PARA_L_TIES = ('AVW_Para_L_fan__PickedSet346_tnof', 'AVW_Para_L_fan__PickedSet347_tnof')
CLUSL = ('CL-L_conn_1', 'CL-L_conn_2', 'CL-L_conn_3', 'CL-R_conn_1', 'CL-R_conn_2', 'CL-R_conn_3',
         'USL-L_conn_1', 'USL-L_conn_2', 'USL-R_conn_1', 'USL-R_conn_2')


def dtol01(m):
    el = m.root.find('Control/solver/dtol')
    old = el.text
    el.text = '0.01'
    m.log.append(f'solver only: dtol {old} -> 0.01 (the displacement-norm tolerance; etol 0.01 still applies)')


def clusl50(m):
    for s in CLUSL:
        spring_scale(m, s, 0.5, CT)


BUILDS = (
    ('L43_lofts_parcus25_avwparatie', 'L38_lofts_parcus25', lambda m: merge_tnof(m, AVW_PARA_L_TIES, loft='AVW-Para-L_fan')),
    ('L43_lofts_parcus25_dtol01', 'L38_lofts_parcus25', dtol01),
    ('L43_springs_parcus10_clusl50', 'L40_springs_parcus10', clusl50),
)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
