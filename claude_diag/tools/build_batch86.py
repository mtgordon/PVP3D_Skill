"""Batch 86 (2026-09-27 ~19:05): the recommended PVW_LA setting (the 30 facets out, penalty 0.5, seg_up 2) checked where it
has not run yet, and the PeB refit on the lofts line (to carry adoptable changes to both lines).
  L82_nopinch_pvwsegup2_pen05       L78_nopinch_pvwsegup2 (source loads) + PVW_LA penalty 5 -> 0.5: the recommended setting
                                    at the source loads; the Ogden-wall control of L80_nopinch_vwyeoh_pvwpen05
  L82_lofts_la3_nopinch_pvwpen05_segup2  L81_lofts_la3_nopinch_pvwpen05 (the lofts line, facets out, penalty 0.5; 3061
                                    iterations) + PVW_LA seg_up 2: does seg_up slow the lofts line at penalty 0.5 too?
  L82_lofts_la3_pebyeoh             L81_lofts_la3_nopinch_pvwpen05 + the Yeoh perineal body (source law refit)
usage: py -3.10 build_batch86.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch77 import pebyeoh

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def pvw(**kw):
    def fn(m):
        assert m.set_contact('PVW_LA', **kw) == 1
    return fn


BUILDS = (('L82_nopinch_pvwsegup2_pen05', 'L78_nopinch_pvwsegup2', pvw(penalty=0.5)),
          ('L82_lofts_la3_nopinch_pvwpen05_segup2', 'L81_lofts_la3_nopinch_pvwpen05', pvw(seg_up=2)),
          ('L82_lofts_la3_pebyeoh', 'L81_lofts_la3_nopinch_pvwpen05', pebyeoh))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
