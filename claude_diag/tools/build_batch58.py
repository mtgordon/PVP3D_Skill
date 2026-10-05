"""Batch 58 (2026-09-26 ~10:30): the candidates and the lofts bases with the source arcus chain mass.

Batch 57: with the source chain mass (x1, not the fast base's x10) the chains move early (mean |u| 4.8 mm at t 0.5 vs 0.6),
the LA ends of the posterior sphincter connectors descend early (-7.6 vs -2.9 mm at t 0.5) and drag the body's posterior
edge down: the springs base turns +0.0 at t 0.5 (fast base +7.3), while P-arcus 10 % + CL/USL 50 % still turns +26.4
(springs) / +27.9 (lofts with the AVW-Para connectors). One change each (the chain density 0.0011 -> 0.00011):
  L51_springs_parcus25_m1                    L38_springs_parcus25
  L51_lofts_parcus25_avwparaconn_m1          L44_lofts_parcus25_avwparaconn
  L51_springs_pc10cu50_sphside50_m1          L48_springs_pc10cu50_sphside50
  L51_lofts_pc10cu50_sphside50_avwparaconn_m1  L49_lofts_pc10cu50_sphside50_avwparaconn
  L51_lofts_la3_pm_m1                        L26_lofts_la3_pm (the lofts reference)
  L51_lofts_avwparaconn_m1                   L46_lofts_avwparaconn
usage: py -3.10 build_batch58.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch57 import mass1

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L51_springs_parcus25_m1', 'L38_springs_parcus25'),
          ('L51_lofts_parcus25_avwparaconn_m1', 'L44_lofts_parcus25_avwparaconn'),
          ('L51_springs_pc10cu50_sphside50_m1', 'L48_springs_pc10cu50_sphside50'),
          ('L51_lofts_pc10cu50_sphside50_avwparaconn_m1', 'L49_lofts_pc10cu50_sphside50_avwparaconn'),
          ('L51_lofts_la3_pm_m1', 'L26_lofts_la3_pm'),
          ('L51_lofts_avwparaconn_m1', 'L46_lofts_avwparaconn'))
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        mass1(m)
        emit(name, m)
