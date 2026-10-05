"""Batch 55 (2026-09-26 ~09:45): does the AVW's mobility turn the body further?

With the AVW-Para lofts as their connectors the lofts line turns the body like the springs line (L46_lofts_avwparaconn +10.1
deg at t 0.66, the lofts base +4.6, the springs base +9.3): the AVW's lateral hold matters to the body's rotation (the
canal contact couples the AVW to the PVW and the body's anterior edge). One change each on L43_springs_parcus10_clusl50
(P-arcus 10 %, CL/USL 50 %; +35.3 deg at t = 1), all NOT IN SOURCE (connective tissue):
  L48_springs_pc10cu50_avwpara50   + AVW-Para-L/R connectors at 50 %
  L48_springs_pc10cu50_pm50        + PM_conn (the 26 perineal-membrane connectors on the distal AVW) at 50 %
  L48_springs_pc10cu50_sphside50   + the side sphincter connectors at 50 %
usage: py -3.10 build_batch55.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CT = 'NOT IN SOURCE (connective tissue)'
BASE = 'L43_springs_parcus10_clusl50'


def avwpara50(m):
    for s in ('AVW-Para-L_conn', 'AVW-Para-R_conn'):
        spring_scale(m, s, 0.5, CT)


BUILDS = (('L48_springs_pc10cu50_avwpara50', avwpara50),
          ('L48_springs_pc10cu50_pm50', lambda m: spring_scale(m, 'PM_conn', 0.5, CT)),
          ('L48_springs_pc10cu50_sphside50', lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.5, CT)))
if __name__ == '__main__':
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, BASE, BASE + '.feb'))
        fn(m)
        emit(name, m)
