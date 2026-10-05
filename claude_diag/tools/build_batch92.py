"""Batch 92 (2026-09-27 ~22:00): the faithful walls' working set carried to the new lofts line.

On the new springs line the faithful (Yeoh) walls with solver rhoi 0 (L85_newline_vwyeoh_rhoi0) got past t 0.74 and
sped up (t 0.771 at 21:56), while without it they crawl at t ~0.66 (L84_springs_newline_vwyeoh, and with the canal's
search_tol 0.1 or knmult 1). The same set on the new lofts line (for faster screens later):
  L87_lofts_newline_vwyeoh_rhoi0   L83_lofts_la3_la100_pebyeoh + the Yeoh walls (all the source's tissue laws) + rhoi 0
and the check that rhoi 0 is a solver setting only (same solution at t = 1 as the line, which runs the Ogden walls):
  L87_springs_newline_rhoi0        L83_springs_la3_la100_pebyeoh + rhoi 0
usage: py -3.10 build_batch92.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch77 import vwyeoh
from build_batch79 import rhoi0

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

BUILDS = (('L87_lofts_newline_vwyeoh_rhoi0', 'L83_lofts_la3_la100_pebyeoh', (vwyeoh, rhoi0)),
          ('L87_springs_newline_rhoi0', 'L83_springs_la3_la100_pebyeoh', (rhoi0,)),
          # 2026-09-28, the user: "go ahead and change rhoi" (adopt rhoi 0 on both lines). The springs line's check above
          # kept its solution and ran 28 % faster; the same check on the lofts line, which becomes its new base:
          ('L87_lofts_newline_rhoi0', 'L83_lofts_la3_la100_pebyeoh', (rhoi0,)))
if __name__ == '__main__':
    for name, base, fns in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        for fn in fns:
            fn(m)
        emit(name, m)
