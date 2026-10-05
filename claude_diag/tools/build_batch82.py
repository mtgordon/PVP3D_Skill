"""Batch 82 (2026-09-27 ~16:45): the faithful materials and the ballooning changes on the lines that now run.

Found since batch 81: the user's option (b) (30 pinch facets out, PVW_LA penalty 5) needs PVW_LA seg_up 2 wherever the
PVW meets the LA at other edge facets (less LA pressure, softer walls): L77_springs_la3_nopinch_pvwsegup2 and
L77_lofts_la3_nopinch_pvwsegup2 reach t = 1 with the old lines' solutions. The Yeoh-wall canal tests of batch 79 lacked it:
at t 0.46 L74_nopinch_vwyeoh crawled both at the canal seam and at the PVW_LA edge onset (node 1984 at 6 m/s).
Ballooning at the source loads (LA area): x2.15 control; x1.60 healthy LA (L74_nopinch_la100); x1.30 at 2x healthy
(L76_nopinch_la200); x1.37 with the arcus chains held along their length (L76_nopinch_chainpin).
Builds (one change each from the named base):
  L78_nopinch_pvwsegup2             L71_aggr_noparcusfacets + PVW_LA seg_up 0 -> 2 (the source-load base with the working
                                    PVW_LA setting; control)
  L78_nopinch_pvwsegup2_vwyeoh      L78_nopinch_pvwsegup2 + the Yeoh vaginal walls (task 1)
  L78_springs_la3_vwyeoh            L77_springs_la3_nopinch_pvwsegup2 (the springs line) + the Yeoh vaginal walls
  L78_springs_la3_pebyeoh           L77_springs_la3_nopinch_pvwsegup2 + the Yeoh perineal body
  L78_lofts_la3_vwyeoh              L77_lofts_la3_nopinch_pvwsegup2 (the lofts line) + the Yeoh vaginal walls (the AVW-Para
                                    loft holds the AVW up: less sliding at the canal seam?)
  L78_springs_la3_la100             L77_springs_la3_nopinch_pvwsegup2 + the healthy LA (NOT IN SOURCE)
  L78_springs_la3_chainpin          L77_springs_la3_nopinch_pvwsegup2 + the arcus chains held (NOT IN SOURCE)
  L78_nopinch_la100_chainpin        L74_nopinch_la100 + the arcus chains held (NOT IN SOURCE; source loads, the combination)
usage: py -3.10 build_batch82.py [NAME ...]   (default: all, in order; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch77 import vwyeoh, pebyeoh
from build_batch78 import la100
from build_batch80 import chainpin

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def pvwsegup2(m):
    assert m.set_contact('PVW_LA', seg_up=2) == 1


BUILDS = (('L78_nopinch_pvwsegup2', 'L71_aggr_noparcusfacets', pvwsegup2),
          ('L78_nopinch_pvwsegup2_vwyeoh', 'L78_nopinch_pvwsegup2', vwyeoh),
          ('L78_springs_la3_vwyeoh', 'L77_springs_la3_nopinch_pvwsegup2', vwyeoh),
          ('L78_springs_la3_pebyeoh', 'L77_springs_la3_nopinch_pvwsegup2', pebyeoh),
          ('L78_lofts_la3_vwyeoh', 'L77_lofts_la3_nopinch_pvwsegup2', vwyeoh),
          ('L78_springs_la3_la100', 'L77_springs_la3_nopinch_pvwsegup2', la100),
          ('L78_springs_la3_chainpin', 'L77_springs_la3_nopinch_pvwsegup2', chainpin),
          ('L78_nopinch_la100_chainpin', 'L74_nopinch_la100', chainpin))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
