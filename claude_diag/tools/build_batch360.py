"""Batch 360 (2026-10-20 night; the user's go (3) of NEXT_SESSION_PROMPT_2026-10-20b: "try allowing PeB to rotate more"). In L352
the perineal body's rotation (about -x) stalls at ~28 deg from t 0.6 while the distal curve span collapses; LA_sphincter_side_conn
(16 springs, scale 0.25 since 10-05) is the only spring family resisting the rotation, its moment growing +3.1 -> +6.5 N mm
(t 0.5 -> 0.8). NOT IN SOURCE (connective tissue, a local spring family, not a wall): the discrete material MAT's <scale> -> SCALE.
usage: py -3.10 build_batch360.py SRC_RUN NEW_RUN [--mat LA_sphincter_side_conn_mat] [--scale 0.1]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new, mat, scale):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    m = re.search(rf'(<discrete_material id="\d+" name="{re.escape(mat)}" type="[^"]+">\s*<scale>)([^<]+)(</scale>)', t)
    old = m.group(2)
    t = t[:m.start()] + m.group(1) + f'{scale:g}' + m.group(3) + t[m.end():]
    msg = f'NOT IN SOURCE: {mat} scale {old} -> {scale:g} (local spring family; tools/build_batch360.py)'
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch360.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: a[a.index(k_) + 1] if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--mat', 'LA_sphincter_side_conn_mat'), float(g('--scale', '0.1')))
