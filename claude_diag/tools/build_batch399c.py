"""Batch 399c (2026-10-22; the user on L399b: the strap wedges tied to the curve, left out to avoid the rest overlap, "are the
exact ones that would prevent the knife edge"). Instead of appending the strap faces to the curve SELF-contact (whose 0.5 mm
search radius reaches through the 0.49 mm strap at the tie), a separate contact: NOT IN SOURCE
  curve_strap: sliding-elastic, primary = curve_all faces (the curve nodes) except faces with a node shared with the strap
  (the tied tip row) or a node behind a strap face within 0.6 mm at rest; secondary = the curve-facing (top) triangles of ALL
  32 strap wedges a side (PM_avw_bottom_left/right_sheet_solid), including the wedges tied to the curve tip, so the surface
  continues across the tied edge; PENALTY PEN auto (as curve_all_self, 5), two_pass 0, search_radius SR 0.25 (< the sheet
  thickness 0.49: the tied tip nodes 0.49 mm behind the far side of their own wedges are out of reach at rest), offset 0.
usage: py -3.10 build_batch399c.py SRC_RUN NEW_RUN [--sr 0.25] [--pen 5]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402

STRAP = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')


def bary_min(p, A, B, C):
    v0, v1, v2 = B - A, C - A, p - A
    d00, d01, d11, d20, d21 = v0 @ v0, v0 @ v1, v1 @ v1, v2 @ v0, v2 @ v1
    den = d00 * d11 - d01 * d01
    v = (d11 * d20 - d01 * d21) / den; w = (d00 * d21 - d01 * d20) / den
    return min(1 - v - w, v, w)


def main(src, new, sr, pen):
    t, X, E = load(src)
    a = t.index('<Surface name="curve_all">'); b = t.index('</Surface>', a)
    F = [[int(v) for v in s.split(',')] for s in re.findall(r'>([0-9,]+)</', t[a:b])]
    strapn = {v for nm, c in E.values() if nm in STRAP for v in c}
    sec = [[c[3], c[4], c[5]] for nm, c in E.values() if nm in STRAP]   # penta6 top = the curve-facing side (build_batch399)
    behind = set()
    for f in sec:
        A, B, C = [X[v] for v in f]; n = np.cross(B - A, C - A); n /= np.linalg.norm(n)
        for v in {v for q in F for v in q} - strapn:
            d = (X[v] - A) @ n
            if -0.6 < d < 0 and bary_min(X[v] - d * n, A, B, C) > -0.05:
                behind.add(v)
    prim = [f for f in F if not set(f) & (strapn | behind)]
    g = lambda f: 'tri3' if len(f) == 3 else 'quad4'   # noqa: E731
    surf = ('\t\t<Surface name="curve_strap_primary">\n' + ''.join(
        f'\t\t\t<{g(f)} id="{k}">{",".join(map(str, f))}</{g(f)}>\n' for k, f in enumerate(prim, 1)) +
        '\t\t</Surface>\n\t\t<Surface name="curve_strap_secondary">\n' + ''.join(
        f'\t\t\t<tri3 id="{k}">{",".join(map(str, f))}</tri3>\n' for k, f in enumerate(sec, 1)) +
        '\t\t</Surface>\n\t\t<SurfacePair name="curve_strap">\n\t\t\t<primary>curve_strap_primary</primary>\n'
        '\t\t\t<secondary>curve_strap_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    a = t.index('<contact name="curve_all_self"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = (t[a:b].replace('name="curve_all_self" surface_pair="curve_all_self"', 'name="curve_strap" surface_pair="curve_strap"'))
    c = re.sub(r'<search_radius>[^<]*</search_radius>', f'<search_radius>{sr:g}</search_radius>', c)
    c = re.sub(r'<penalty>[^<]*</penalty>', f'<penalty>{pen:g}</penalty>', c)
    for s_ in ('<two_pass>0</two_pass>', '<offset>0</offset>', '<laugon>PENALTY</laugon>'):
        assert s_ in c, s_
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE (the user, 10-22): contact curve_strap (sliding-elastic, as curve_all_self but search_radius {sr:g}, '
           f'penalty {pen:g} auto, two_pass 0, offset 0): primary = {len(prim)} of {len(F)} curve_all faces (left out: faces with a '
           f'node shared with the strap, or with one of {len(behind)} curve nodes behind a strap face within 0.6 mm at rest); '
           f'secondary = the curve-facing top triangles of all {len(sec)} strap wedges (both sides, the wedges tied to the curve '
           f'tip included): the curve nodes meet the strap surface right from the tied edge; curve_all_self unchanged')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch399c.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--sr', 0.25), g('--pen', 5))
