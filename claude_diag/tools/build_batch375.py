"""Batch 375 (2026-10-21; against the lift-and-slam of L359: the curve pressure squeezes the whole soft curve, which lengthens
and jacks the cervix ~22 mm up, then drops 42 mm at t ~0.9). The hairpin it cures is only in the distal curve, so: the curve
pressure on the DISTAL faces only. NOT IN SOURCE (the source has no curves; the curve pressure is ours): a new surface
Load-curves_distal_surf = the Load-curves_surf faces whose centre is at rest z < ZCUT, and the Load-curves surface_load uses it
(same 0.014 MPa, lc 1). The full Load-curves_surf stays defined (unused).
usage: py -3.10 build_batch375.py SRC_RUN NEW_RUN [--zcut -12]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new, zcut):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    ms = re.search(r'\t\t<Surface name="Load-curves_surf">\n(.*?)\t\t</Surface>\n', t, re.S)
    rows = re.findall(r'\t\t\t<(\w+) id="\d+">([0-9,]+)</\1>\n', ms.group(1))
    keep = [(ty, c) for ty, c in rows if np.mean([X[int(v)][2] for v in c.split(',')]) < zcut]
    blk = ('\t\t<Surface name="Load-curves_distal_surf">\n' + ''.join(f'\t\t\t<{ty} id="{i}">{c}</{ty}>\n' for i, (ty, c) in enumerate(keep, 1))
           + '\t\t</Surface>\n')
    t = t[:ms.end()] + blk + t[ms.end():]
    old = '<surface_load name="Load-curves" surface="Load-curves_surf" type="pressure">'
    assert t.count(old) == 1
    t = t.replace(old, '<surface_load name="Load-curves" surface="Load-curves_distal_surf" type="pressure">')
    msg = (f'NOT IN SOURCE (the curve pressure is ours): the curve pressure Load-curves on the DISTAL faces only: '
           f'Load-curves_distal_surf = {len(keep)} of the {len(rows)} Load-curves_surf faces with centre at rest z < {zcut:g} '
           f'(0.014 MPa, lc 1 unchanged) (tools/build_batch375.py)')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch375.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--zcut') + 1]) if '--zcut' in a else -12.0)
