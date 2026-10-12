"""Batch 191 (2026-10-06, the user: "yes, start option 1"). L189 (curve pressure, L169 base) flapped from t ~0.68: the
curves' nodes (mid-length, rest height -5..+11 mm) sped up ~10x per 0.02 t, the curves' AVW half most, until the run
crawled (p95 node speed 0.06 -> 2-4 mm per 0.001 t; L169 without the curve pressure stayed at 0.04-0.07). One change on
L190_botfull_curvepress_P1 (the user's equal PM-arch springs + the curve pressure):
  L191_botfull_curvepress_cdamp1000_P1   the curves' mass damping (wall_damping_canal_wrap_hex_avw / _pvw) C 200 -> 1000
                                         (NOT IN SOURCE; the walls keep C 200)
usage: py -3.10 build_batch191.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L190_botfull_curvepress_P1'
NAME = 'L191_botfull_curvepress_cdamp1000_P1'

if __name__ == '__main__':
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    for d in ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw'):
        a = t.index(f'<body_load name="wall_damping_{d}"')
        b = t.index('</body_load>', a)
        blk, k = re.subn(r'<C>200</C>', '<C>1000</C>', t[a:b])
        assert k == 1, d
        t = t[:a] + blk + t[b:]
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch191.py)\n'
                   f'NOT IN SOURCE: the curves\' mass damping (canal_wrap_hex_avw / _pvw) C 200 -> 1000 /s (walls keep 200)\n')
