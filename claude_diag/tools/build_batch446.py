"""Batch 446 (2026-10-23; the user: "build and start L446 = L444 + the strap material x1/3"). NOT IN SOURCE (our strap): the fan
strap materials PM_avw_bottom_left / right_sheet_mat get their matrix c1 and fibre E scaled by FACTOR (k unchanged), as L431
did to L423b by hand (c1 0.0833 -> 0.0277667, E 29.5 -> 9.83333); L431 removed L423b's late hold collapse.
usage: py -3.10 build_batch446.py SRC_RUN NEW_RUN [--factor 0.333333]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from build_batch162 import write  # noqa: E402

MATS = ('PM_avw_bottom_left_sheet_mat', 'PM_avw_bottom_right_sheet_mat')


def main(src, new, f):
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8').read()
    rep = []
    for name in MATS:
        m = re.search(r'<material id="\d+" name="%s"[^>]*>.*?</material>' % name, t, re.S)
        assert m, name
        old = {}

        def sc(mm):
            v = float(mm.group(2)); old.setdefault(mm.group(1), (v, v * f))
            return f'<{mm.group(1)}>{v * f:g}</{mm.group(1)}>'
        blk = re.sub(r'<(c1|E)>([^<]+)</\1>', sc, m.group(0))
        t = t[:m.start()] + blk + t[m.end():]
        rep.append(f'{name}: ' + ', '.join(f'{k} {a:g} -> {b:g}' for k, (a, b) in old.items()))
    msg = f'NOT IN SOURCE (our strap): the fan strap material x{f:g} (as L431; k unchanged): ' + '; '.join(rep) + '. One change.'
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch446.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--factor') + 1]) if '--factor' in a else 1 / 3)
