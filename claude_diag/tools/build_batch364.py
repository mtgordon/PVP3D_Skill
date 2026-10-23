"""Batch 364 (2026-10-21 night; NEXT_SESSION_PROMPT_2026-10-20b: carry what helps against the curve collapse to the full candidate).
L359 (L352 + the source's curve pressure) kept the distal curve stretched (span 1.2x rest) and its self-penetration <= 0.46 mm
through t 0.84 where L352 snapped to 0.8 mm at 0.67-0.70 and 0.9 mm at 0.83. NOT IN SOURCE (corrected 10-21: the source has no curves; our wall pressure carried onto them): Load-curves on
Load-curves_surf (512 outer curve faces, 0.014 MPa, lc 1, as Load-AVW), removed on the tube line in L258, added after the last
surface_load of SRC (Load-curves_surf must already be in SRC's mesh).
usage: py -3.10 build_batch364.py SRC_RUN NEW_RUN"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
BLOCK = ('\t\t<surface_load name="Load-curves" surface="Load-curves_surf" type="pressure">\n'
         '\t\t\t<pressure lc="1">0.014</pressure>\n\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n'
         '\t\t\t<linear>0</linear>\n\t\t\t<shell_bottom>0</shell_bottom>\n\t\t</surface_load>\n')


def main(src, new):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    assert '<Surface name="Load-curves_surf">' in t and 'name="Load-curves"' not in t
    i = t.rindex('</surface_load>\n') + len('</surface_load>\n')
    t = t[:i] + BLOCK + t[i:]
    msg = ('NOT IN SOURCE (the source has no curves; our wall pressure carried onto them): the curve pressure Load-curves (Load-curves_surf, 512 outer curve faces, 0.014 MPa, lc 1, as '
           'Load-AVW) back on (removed on the tube line, L258); as L359 on L352 (tools/build_batch364.py)')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch364.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
