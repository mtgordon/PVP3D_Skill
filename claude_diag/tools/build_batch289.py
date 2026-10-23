"""Batch 289 (2026-10-15; the user: "yes, start the damping run").
The curves' fold-through at MAX is a fast collapse (L287_dtwin_off06_MAX: steps ~0.00013 through it, yet curve nodes move
0.3-0.7 mm per step). Here, to slow it physically (NOT IN SOURCE; a damping change only): every curve hex
(canal_wrap_hex_avw / _pvw) with its centre at rest z < ZCUT moves to its own domain *_distal (same material Curve_soft), and the
distal domains + canal_wrap_hex_pvw_fine (all of it distal) get mass damping C = CD instead of 1000; the rest unchanged.
usage: py -3.10 build_batch289.py SRC_RUN NEW_RUN [--c 10000] [--zcut -15]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, cd=10000.0, zcut=-15.0):
    t, X, E = load(src)
    moved = {}
    for dom in ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw'):
        ids = [k for k, (n, c) in E.items() if n == dom and np.mean([X[v][2] for v in c]) < zcut]
        a = t.index(f'<Elements type="hex8" name="{dom}">'); b = t.index('</Elements>', a)
        blk = t[a:b]
        rows = []
        for k in ids:
            row = f'\t\t\t<elem id="{k}">' + ','.join(map(str, E[k][1])) + '</elem>\n'
            assert row in blk, k
            blk = blk.replace(row, '', 1); rows.append(row)
        t = t[:a] + blk + t[b:]
        moved[dom] = rows
    a = t.index('<Elements type="hex8" name="canal_wrap_hex_pvw_fine">'); b = t.index('</Elements>', a) + len('</Elements>\n')
    t = t[:b] + ''.join(f'\t\t<Elements type="hex8" name="{d}_distal">\n' + ''.join(r) + '\t\t</Elements>\n'
                        for d, r in moved.items()) + t[b:]
    for d in moved:
        t = t.replace(f'\t\t<SolidDomain name="{d}" mat="Curve_soft" />\n',
                      f'\t\t<SolidDomain name="{d}" mat="Curve_soft" />\n\t\t<SolidDomain name="{d}_distal" mat="Curve_soft" />\n', 1)
    m = re.search(r'\t\t<body_load name="wall_damping_canal_wrap_hex_pvw_fine" .*?</body_load>\n', t, re.S)
    fine = m.group(0)
    assert '<C>1000</C>' in fine
    new_bl = fine.replace('<C>1000</C>', f'<C>{cd:g}</C>') + ''.join(
        fine.replace('canal_wrap_hex_pvw_fine', f'{d}_distal').replace('<C>1000</C>', f'<C>{cd:g}</C>') for d in moved)
    t = t.replace(fine, new_bl, 1)
    n = {d: len(r) for d, r in moved.items()}
    msg = (f'NOT IN SOURCE (damping only): curve hexes with centre at rest z < {zcut} moved to their own domains {n} '
           f'(*_distal, Curve_soft); those + canal_wrap_hex_pvw_fine get mass damping C {cd:g} (was 1000); the rest unchanged')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch289.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--c', 10000.0), g('--zcut', -15.0))
