"""Batch 49 (2026-09-26): the lofts line with the distal P-arcus cut, plus the next hold (the side sphincters).

Batch 46: on the lofts line every variant that weakens the proximal P-arcus crawls at t ~0.43-0.52 (parcus50 got through
in 63 min; parcus25, parcusprox0 and the USL-L-merged parcus50 crawl there), but cutting only the 12 distal connectors
runs smoothly (L38_lofts_parcusdist0: t 0.76 at 34 min, 1 failed, +16.1 deg). One change:
  L41_lofts_parcusdist0_sphside50   L38_lofts_parcusdist0 + side sphincters at 50 %    NOT IN SOURCE (connective tissue)
  L41_springs_parcusdist0_sphside50 L38_springs_parcusdist0 + the same (08:33; its springs twin)
  L42_lofts_parcusdist0_nosphside   L38_lofts_parcusdist0 + side sphincters cut          DIAGNOSTIC (08:35)
usage: py -3.10 build_batch49.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CT = 'NOT IN SOURCE (connective tissue)'
BUILDS = (
    ('L41_lofts_parcusdist0_sphside50', 'L38_lofts_parcusdist0', lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.5, CT)),
    ('L41_springs_parcusdist0_sphside50', 'L38_springs_parcusdist0',
     lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.5, CT)),
    ('L42_lofts_parcusdist0_nosphside', 'L38_lofts_parcusdist0', lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.0)),
)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
