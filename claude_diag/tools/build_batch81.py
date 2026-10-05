"""Batch 81 (2026-09-27 ~09:05): the user's option (b) on the springs line (Load-LA 1/3) crawls at the PVW_LA onset.

L74_springs_la3_nopinch (L26_springs_la3_pm + the 30 pinch facets out + PVW_LA penalty 0.5 -> 5) matched
L26_springs_la3_pm exactly to t 0.50, then crawled where the PVW first touches the LA: on the PVW contact surface's
lateral edge row (nodes 1898/1899; not P-arcus connector ends), node 1898 reversing direction. L74_nopinch_la0 (no LA
pressure) stalled the same way at t 0.435 (edge nodes 1983/1984, 1899/5402, 7327). With less pressure on the LA it sits
higher and the PVW meets it earlier, at edge facets other than the 30 removed. At penalty 0.5 (L26_springs_la3_pm) the
same onset passes. One change each from L74_springs_la3_nopinch:
  L77_springs_la3_nopinch_pvwsegup2  PVW_LA seg_up 0 -> 2 (Test B: seg_up 2 got past the first onset in L70_aggr_segup2,
                                     then crawled at the P-arcus pinch facets, which option (b) removes)
  L77_springs_la3_nopinch_pvwpen05   PVW_LA penalty 5 -> 0.5 (= L26_springs_la3_pm with the 30 facets out)
usage: py -3.10 build_batch81.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def pvw(**kw):
    def fn(m):
        assert m.set_contact('PVW_LA', **kw) == 1
    return fn


BUILDS = (("L77_springs_la3_nopinch_pvwsegup2", pvw(seg_up=2)),
          ('L77_springs_la3_nopinch_pvwpen05', pvw(penalty=0.5)))
if __name__ == '__main__':
    base = 'L74_springs_la3_nopinch'
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
# Added 09:56 (built inline with the same operation): L77_lofts_la3_nopinch_pvwsegup2 = L74_lofts_la3_nopinch (the lofts
# line with option (b), which crawled at the same PVW edge onset at t 0.467) + PVW_LA seg_up 0 -> 2.
