"""Batch 124 (2026-09-29, the user's choice "#2": adopt the rounded tube L123_tube_round_r30 and run the paper's cases on it).
The paper's cases P1 / P2 (Luo et al. 2015; build_batch107.py's P1 / P2: levator 80 / 40 % of healthy, apical CL / USL
x 0.7 / 0.4, Parcus x 0.15) on the tube and, like for like, on the springs line's base it was built from.
  L124_tube_round_r30_paperP1 / _paperP2            L123_tube_round_r30 + P1 / P2
  L124_newline_vwyeoh_rhoi0_seamspr001_paperP1 / P2 L91_newline_vwyeoh_rhoi0_seamspr001 + P1 / P2
usage: py -3.10 build_batch124.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch107 import P1, P2  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
TUBE, BASE = 'L123_tube_round_r30', 'L91_newline_vwyeoh_rhoi0_seamspr001'
BUILDS = (('L124_tube_round_r30_paperP1', TUBE, P1), ('L124_tube_round_r30_paperP2', TUBE, P2),
          ('L124_newline_vwyeoh_rhoi0_seamspr001_paperP1', BASE, P1),
          ('L124_newline_vwyeoh_rhoi0_seamspr001_paperP2', BASE, P2))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
