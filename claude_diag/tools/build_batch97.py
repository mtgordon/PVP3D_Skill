"""Batch 97 (2026-09-28 ~10:50): does closing the canal seam bring the anterior wall toward the paper? (the user's request)

With the seam free, FEBio's anterior wall sits 13-20 mm lower than the paper's published cases (Luo et al. 2015: Ba above
the hymen in every case; L68_pen5_paperP1: Ba +13.3 mm, Bp -0.8 vs the paper's +4). On the springs line a seam spring of
0.01 N/mm per node (NOT IN SOURCE) lifts the anterior wall 11 mm (L91_springs_newline_rhoi0_seamspr001: Ba -14.6 vs -3.3)
and is the softest seam connection that lets the faithful walls converge. One change:
  L92_pen5_paperP1_seamspr001   L68_pen5_paperP1 (the paper's case P1 on L19_pen5, the Abaqus replica: levator 20 %
                                impaired, apical 30 %, posterior 85 %; source loads, Ogden walls) + the 0.01 N/mm seam springs
usage: py -3.10 build_batch97.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch90 import seam_tie

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

BUILDS = (('L92_pen5_paperP1_seamspr001', 'L68_pen5_paperP1', lambda m: seam_tie(m, penalty=0.01, maxaug=0)),
          # 11:20, the user: "run P2 with the springs too" (P1: Ba +13.3 -> -4.6, the paper's side; Bp -0.8 -> -8.7,
          # further from the paper's +4). L68_pen5_paperP2: levator 60 % impaired, apical 60 %: Ba +20.3, Bp +3.6 (paper +9)
          ('L92_pen5_paperP2_seamspr001', 'L68_pen5_paperP2', lambda m: seam_tie(m, penalty=0.01, maxaug=0)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        print(name, ':', m.log[-1][:120])
        emit(name, m)
