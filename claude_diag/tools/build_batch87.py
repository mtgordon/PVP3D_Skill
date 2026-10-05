"""Batch 87 (2026-09-27 ~20:05): the new lines with the user's decisions of ~20:00.

The user: (1) PVW_LA on both lines = the 30 pinch facets out + penalty 0.5 + seg_up 2 (the bases
L81_springs_la3_nopinch_pvwsegup2_pen05 and L82_lofts_la3_nopinch_pvwpen05_segup2, both t = 1 with the old lines'
solutions); (3) the healthy LA from the source's own PCM-LA_Yamada100% table AND Load-LA 1/3 (already on the lines);
adopt the perineal-body refit (the source's Marlow data fit to Yeoh); (4) keep pushing the faithful walls. One change each:
  L83_springs_la3_la100            L81_springs_la3_nopinch_pvwsegup2_pen05 + the healthy LA (NOT IN SOURCE: the source uses
                                   the 50 % law; the 100 % table is defined in the .inp, unused)
  L83_lofts_la3_la100              L82_lofts_la3_nopinch_pvwpen05_segup2 + the healthy LA
  L83_springs_la3_la100_pebyeoh    L83_springs_la3_la100 + the Yeoh perineal body: THE NEW SPRINGS LINE
  L83_lofts_la3_la100_pebyeoh      L83_lofts_la3_la100 + the Yeoh perineal body: THE NEW LOFTS LINE
  L84_springs_newline_vwyeoh       L83_springs_la3_la100_pebyeoh + the Yeoh vaginal walls (all the source's tissue laws;
                                   the reference for the canal-contact work)
usage: py -3.10 build_batch87.py [NAME ...]   (default: all, in order; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch77 import vwyeoh, pebyeoh
from build_batch78 import la100

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

BUILDS = (('L83_springs_la3_la100', 'L81_springs_la3_nopinch_pvwsegup2_pen05', la100),
          ('L83_lofts_la3_la100', 'L82_lofts_la3_nopinch_pvwpen05_segup2', la100),
          ('L83_springs_la3_la100_pebyeoh', 'L83_springs_la3_la100', pebyeoh),
          ('L83_lofts_la3_la100_pebyeoh', 'L83_lofts_la3_la100', pebyeoh),
          ('L84_springs_newline_vwyeoh', 'L83_springs_la3_la100_pebyeoh', vwyeoh))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
