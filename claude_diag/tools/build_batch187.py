"""Batch 187 (2026-10-06, the user, on the tube's PM_avw_bottom springs being ~1e5 x softer than PM_conn_mat: "make them all
the stiffer value").
In the tube runs each curve's distal end has the source's 2 PM_avw_bottom springs a side (law PM-LA x1e-5, i.e. off) and the
3 added curve_end_PM_conn springs a side (PM_conn_mat, full; NOT IN SOURCE, L166). One change on the working unrefined tube run
L169_wrapfull_endspr_selfall_walls035_P1 (self-contact penalty 5; fold 2 / 22 / 19 / 26 nodes inside at t 0.8 / 0.85 / 0.9 / 1):
  L187_selfall_botfull_P1   PM_avw_bottom_left_conn_mat and PM_avw_bottom_right_conn_mat get PM_conn_mat's force points
                            (x1e5; NOT IN SOURCE / NOT FAITHFUL: the source has them off)
Judge with tools/curve_self_pen.py RUN 0.8,0.85,0.9,1 and summary_table.py.
usage: py -3.10 build_batch187.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L169_wrapfull_endspr_selfall_walls035_P1'


def block(t, name):
    a = t.index(f'name="{name}" type="nonlinear spring">')
    return a, t.index('</discrete_material>', a)


if __name__ == '__main__':
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    a, b = block(t, 'PM_conn_mat')
    pts = re.search(r'<points>.*?</points>', t[a:b], re.S).group(0)
    for side in ('left', 'right'):
        a, b = block(t, f'PM_avw_bottom_{side}_conn_mat')
        new, k = re.subn(r'<points>.*?</points>', pts, t[a:b], flags=re.S)
        assert k == 1
        t = t[:a] + new + t[b:]
    write('L187_selfall_botfull_P1', t,
          f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch187.py)\n'
          f'NOT IN SOURCE / NOT FAITHFUL (the user, 2026-10-06): PM_avw_bottom_left_conn_mat and PM_avw_bottom_right_conn_mat '
          f'force points replaced by PM_conn_mat\'s (source law PM-LA x1e-5 -> PM_conn_mat, x~1e5; 4 springs, 2 a side, on the '
          f'curves\' distal ends)\n')
