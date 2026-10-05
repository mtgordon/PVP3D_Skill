"""Batch 108 (2026-09-28 ~22:00, the user's choice): the distal anterior wall's support. NOT IN SOURCE as given.

The canal seam's most-open node (1773, the AVW's distal side corner at the introitus) is also one end of the PM_avw_bottom
connectors (2 per side: AVW nodes 1773 / 1821 on the left, to the perineal membrane). In the source .inp their table is
1/100,000 of the PM connector curve (6e-05 N at 4.7 mm, 0.003 N at 47 mm vs PM 5.8 / 296 N): the distal AVW corners are
free, where the anterior wall bulges out through the introitus and peels off the PVW, and where FEBio's Ba sits 13-32 mm
below the hymen in the paper's cases (the paper: above it). Which impairment baseline the paper used is unknown (the user
does not know which version of the .inp they were given). Test: PM_avw_bottom_left / right_conn scale 1 -> 1e5 (= the PM
curve, "normal" support if PM is the baseline), one change each from:
  L108_springs_newline_rhoi0_avwbot1e5             L87_springs_newline_rhoi0 (the springs line, the standard case)
  L108_springs_newline_rhoi0_paperP1_avwbot1e5     L107_springs_newline_rhoi0_paperP1 (the paper's P1, Ogden walls)
  L108_springs_newline_rhoi0_paperP2_avwbot1e5     L107_springs_newline_rhoi0_paperP2 (the paper's P2)
  L108_newline_vwyeoh_rhoi0_avwbot1e5              L85_newline_vwyeoh_rhoi0 (the faithful walls; crawl at t 0.78)
Bracket: the same three springs-line cases at x1e4 (the AVW-Para level), *_avwbot1e4.
usage: py -3.10 build_batch108.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch44 import spring_scale

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
LABEL = 'NOT IN SOURCE (distal anterior support: PM_avw_bottom at the PM curve\'s strength)'


def avwbot(s):
    def fn(m):
        for n in ('PM_avw_bottom_left_conn', 'PM_avw_bottom_right_conn'):
            spring_scale(m, n, s, LABEL)
    return fn


BUILDS = (('L108_springs_newline_rhoi0_avwbot1e5', 'L87_springs_newline_rhoi0', avwbot(1e5)),
          ('L108_springs_newline_rhoi0_paperP1_avwbot1e5', 'L107_springs_newline_rhoi0_paperP1', avwbot(1e5)),
          ('L108_springs_newline_rhoi0_paperP2_avwbot1e5', 'L107_springs_newline_rhoi0_paperP2', avwbot(1e5)),
          ('L108_newline_vwyeoh_rhoi0_avwbot1e5', 'L85_newline_vwyeoh_rhoi0', avwbot(1e5)),
          # the bracket (2026-09-28 night, NEXT_SESSION_PROMPT_2026-09-28.md: "you may bracket with x1e4 (the AVW-Para
          # level)"): the same three cases at 1/10 of the PM curve
          ('L108_springs_newline_rhoi0_avwbot1e4', 'L87_springs_newline_rhoi0', avwbot(1e4)),
          ('L108_springs_newline_rhoi0_paperP1_avwbot1e4', 'L107_springs_newline_rhoi0_paperP1', avwbot(1e4)),
          ('L108_springs_newline_rhoi0_paperP2_avwbot1e4', 'L107_springs_newline_rhoi0_paperP2', avwbot(1e4)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
