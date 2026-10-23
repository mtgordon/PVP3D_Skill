"""Batch 399 (2026-10-22; the user: "add the surface of the pm_avw strap to the contact list. That will eliminate the edge of the
element because it's attached to the strap"). In L394 the curve penetrates at the distal end corner where the strap is tied
straight onto the curve tip (face 25892 / 25893 / 25898 / 25897 of canal_wrap_hex_pvw_fine 28676: two nodes shared with the strap,
in curve_all; the thin end hexes there swell 2.5x / crush 0.3x and curve nodes from higher up slip past the tied edge).
NOT IN SOURCE: the strap's (PM_avw_bottom_left/right_sheet_solid, 32 wedges a side, one layer ~0.49 mm thick) exterior
triangle faces on the side facing the curve are appended to curve_all (the curve self-contact, curve_all_self), so the
self-contact surface continues past the tied edge onto the strap. Only ONE side of the strap: the sheet is thinner than
curve_all_self's search_radius 0.5 (gotcha 30: a self-contact reaching across the part's own thickness grabs its opposite face
at rest). The side is the one whose outward normal points toward the curve nodes within DIST mm (median over the sheet).
The strap's edge quads are left out.
usage: py -3.10 build_batch399.py SRC_RUN NEW_RUN [--dist 15] [--all (keep the wedges tied to the curve)] [--sr R (curve_all_self search_radius)]"""
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402

STRAP = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
TRI = ((0, 2, 1), (3, 4, 5))   # penta6 bottom / top, outward


def main(src, new, dist, keep_all=False, sr=None):
    t, X, E = load(src)
    cn = np.array(sorted({v for nm, c in E.values() if nm in CURVES for v in c}))
    tree = cKDTree(np.array([X[v] for v in cn]))
    add, notes = [], []
    curvn = {v for nm, c in E.values() if nm in CURVES for v in c}
    for s in STRAP:
        wed = [c for nm, c in E.values() if nm == s]
        assert wed and all(len(c) == 6 for c in wed), s
        # (L399b) the wedges sharing nodes with the curves sit in the curve's surface at rest: their faces / nodes lie
        # 0.49 mm behind curve facets (L399 kicked at the first step); left out, the strap facets start ~1.8 mm from the tie
        nall = len(wed)
        if not keep_all:
            wed = [c for c in wed if not set(c) & curvn]
        score = []
        for side in (0, 1):
            d = []
            for c in wed:
                q = [c[i] for i in TRI[side]]; P = np.array([X[v] for v in q]); ctr = P.mean(0)
                n = np.cross(P[1] - P[0], P[2] - P[0]); n /= np.linalg.norm(n)
                near = tree.query_ball_point(ctr, dist)
                if near:
                    d.append(np.median((np.array([X[cn[j]] for j in near]) - ctr) @ n))
            score.append(np.median(d) if d else -1e9)
        side = int(np.argmax(score))
        faces = [[c[i] for i in TRI[side]] for c in wed]
        add += faces
        notes.append(f'{s}: {len(faces)} {"top" if side else "bottom"} triangles ({nall - len(wed)} of {nall} wedges sharing '
                     f'nodes with the curves left out; median curve-ward distance along the '
                     f'outward normal {score[side]:+.1f} mm vs {score[1 - side]:+.1f} for the other side)')
    a = t.index('<Surface name="curve_all">'); b = t.index('\t\t</Surface>', a)
    n0 = max(int(v) for v in re.findall(r'id="(\d+)"', t[a:b]))
    old = {frozenset(int(v) for v in s.split(',')) for s in re.findall(r'>([0-9,]+)</', t[a:b])}
    assert not any(frozenset(f) in old for f in add)
    blk = ''.join(f'\t\t\t<tri3 id="{n0 + k}">{",".join(map(str, f))}</tri3>\n' for k, f in enumerate(add, 1))
    t = t[:b] + blk + t[b:]
    if sr is not None:   # (L399d, the user) a smaller curve_all_self search radius, below the 0.49 mm sheet
        a = t.index('<contact name="curve_all_self"'); b = t.index('</contact>', a)
        c = re.sub(r'<search_radius>[^<]*</search_radius>', f'<search_radius>{sr:g}</search_radius>', t[a:b])
        assert c != t[a:b]
        t = t[:a] + c + t[b:]
    msg = (f'NOT IN SOURCE (the user, 10-22): the strap faces facing the curve appended to curve_all (the curve self-contact '
           f'curve_all_self; facets {n0 + 1}-{n0 + len(add)}): ' + '; '.join(notes) +
           '; edge quads and the far side left out (sheet ~0.49 mm)' + (f'; curve_all_self search_radius 0.5 -> {sr:g}' if sr is not None else '; contact settings unchanged'))
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch399.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--dist') + 1]) if '--dist' in a else 15.0, '--all' in a,
         float(a[a.index('--sr') + 1]) if '--sr' in a else None)
