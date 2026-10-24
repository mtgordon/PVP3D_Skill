"""Batch 453 (2026-10-23; the user: the strap should swing out with the AVW like the PM connectors to the AVW do). NOT IN SOURCE
(our strap): the strap weighs ~75 mg (inertia / mass damping negligible); what stops it swinging is its matrix (Mooney-Rivlin c1
0.0278 after x1/3, ~13x the wall's) and k 10 (~29x the wall's), which resist shear and bending. Here the fan strap materials
(PM_avw_bottom_left / right_sheet_mat, also on the _wide strips) keep their tension-only fibres and get the matrix c1 x FACTOR
and k = K, so the strap carries tension along its fibres like a sheet of connectors and has almost no shear / bending stiffness.
usage: py -3.10 build_batch453.py SRC_RUN NEW_RUN [--c1factor 0.01] [--k 0.1]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from build_batch162 import write  # noqa: E402

MATS = ('PM_avw_bottom_left_sheet_mat', 'PM_avw_bottom_right_sheet_mat')


def main(src, new, f, k):
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8').read()
    rep = []
    for name in MATS:
        m = re.search(r'<material id="\d+" name="%s"[^>]*>.*?</material>' % name, t, re.S)
        assert m, name
        blk = m.group(0)
        c1 = float(re.search(r'<c1>([^<]+)</c1>', blk).group(1)); k0 = float(re.search(r'<k>([^<]+)</k>', blk).group(1))
        blk = re.sub(r'<c1>[^<]+</c1>', f'<c1>{c1 * f:g}</c1>', blk, count=1)
        blk = re.sub(r'<k>[^<]+</k>', f'<k>{k:g}</k>', blk, count=1)
        t = t[:m.start()] + blk + t[m.end():]
        rep.append(f'{name}: matrix c1 {c1:g} -> {c1 * f:g}, k {k0:g} -> {k:g}')
    msg = ('NOT IN SOURCE (our strap): connector-like strap (fibres unchanged, matrix x%g, k %g: almost no shear / bending '
           'stiffness, so it swings out like the PM connectors): ' % (f, k) + '; '.join(rep) + '. One change.')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch453.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    opt = lambda n, d: (a[a.index(n) + 1] if n in a else d)  # noqa: E731
    main(a[0], a[1], float(opt('--c1factor', 0.01)), float(opt('--k', 0.1)))
