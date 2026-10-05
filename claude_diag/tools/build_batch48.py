"""Batch 48 (2026-09-26): past P-arcus, the side sphincter connectors are the next hold on the body's rotation.

With all P-arcus connectors cut (L36_springs_noparcus) the rotation levels off at ~45 deg from t ~0.55: the anterior edge
stops at -23 mm (t 0.82) and the side sphincter connectors, which carried little in the base, now resist (-31 N mm at
t = 1, 8.4 N up, up to 9 mm stretch). Batch 46 (springs, t 0.66): P-arcus 25 % +24.3 deg, distal half cut +21.3, proximal
half cut +18.9, all cut +46.1 (base +9.3): the P-arcus hold is spread along the PVW. One change each, springs line:
  L40_springs_parcus25_sphside50     L38_springs_parcus25 + side sphincters at 50 %      NOT IN SOURCE (connective tissue)
  L40_springs_parcus10               L26_springs_la3_pm + P-arcus at 10 %                 NOT IN SOURCE (connective tissue)
  L40_springs_noparcus_sphside50     L36_springs_noparcus + side sphincters at 50 %       (on a DIAGNOSTIC base)
  L40_springs_noparcus_nosphside     L36_springs_noparcus + side sphincters cut           DIAGNOSTIC
usage: py -3.10 build_batch48.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CT = 'NOT IN SOURCE (connective tissue)'
BUILDS = (
    ('L40_springs_parcus25_sphside50', 'L38_springs_parcus25', lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.5, CT)),
    ('L40_springs_parcus10', 'L26_springs_la3_pm', lambda m: spring_scale(m, 'Parcus_conn', 0.1, CT)),
    ('L40_springs_noparcus_sphside50', 'L36_springs_noparcus',
     lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.5, CT)),
    ('L40_springs_noparcus_nosphside', 'L36_springs_noparcus',
     lambda m: spring_scale(m, 'LA_sphincter_side_conn', 0.0)),
)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
