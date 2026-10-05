"""Batch 84 (2026-09-27 ~18:05): carry PVW_LA penalty 0.5, the best fix for the Yeoh walls so far.

At 18:03 L79_springs_la3_vwyeoh_pvwpen05 (the springs line + the Yeoh walls + PVW_LA penalty 5 -> 0.5) went t 0.22 -> 0.61
in 30 min with 4 failed attempts, against the control L78_springs_la3_vwyeoh 0.43 -> 0.50; the canal penalty 0.5
(L79_springs_la3_vwyeoh_se1pen05) and the 26 PVW edge facets out (L79_springs_la3_vwyeoh_pvwedge6) were in between.
  L80_springs_la3_vwyeoh_pvwpen05_pebyeoh  L79_springs_la3_vwyeoh_pvwpen05 + the Yeoh perineal body: the springs line
                                           with all the source's tissue laws (walls, PeB, LA) refit
  L80_nopinch_vwyeoh_pvwpen05              L78_nopinch_pvwsegup2_vwyeoh (source loads; crawled past t 0.6 at +0.03 per
                                           30 min) + PVW_LA penalty 5 -> 0.5
usage: py -3.10 build_batch84.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch77 import pebyeoh

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def pvwpen05(m):
    assert m.set_contact('PVW_LA', penalty=0.5) == 1


BUILDS = (('L80_springs_la3_vwyeoh_pvwpen05_pebyeoh', 'L79_springs_la3_vwyeoh_pvwpen05', pebyeoh),
          ('L80_nopinch_vwyeoh_pvwpen05', 'L78_nopinch_pvwsegup2_vwyeoh', pvwpen05))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
