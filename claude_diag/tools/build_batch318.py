"""Batch 318 (2026-10-17 night; the loft-contacts plan step 2 (memory project-loft-contacts-plan) for the SOLID lofts of
build_batch312: AVW-Para <-> every vaginal wall part. L297 (shell sheets) could not hold: one-sided shell contact and the
curve-end faces turning toward the sheet; a solid loft has outward faces on every side).
NOT IN SOURCE: sliding-elastic contact AVW_Para_solid_walls:
  primary: the outward faces (top, bottom, free edges) of the AVW-Para solid-loft wedges (AVW-Para-L/R_sheet_solid) that have
           no node on the wall surface, minus wedges within one wedge of the attached row (--ring 1) and any face with a node
           inside a wall hex at rest (reported)
  secondary: every wall boundary face (build_batch297 WALLS) within DMAX mm of a primary node at rest, except the lumen faces
           (SlidingElastic1Primary / Secondary)
Settings as walls_LA but penalty PEN, two_pass 0, search_radius SR (0.5: below the wall thickness; gotcha 30), offset 0.
usage: py -3.10 build_batch318.py SRC_RUN NEW_RUN [--sr 0.5] [--pen 0.5] [--dmax 10] [--ring 1]"""
import os
import re
import sys
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch297 import WALLS, HF  # noqa: E402

SOL = ('AVW-Para-L_sheet_solid', 'AVW-Para-R_sheet_solid')
T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]


def main(src, new, sr, pen, dmax, ring):
    t, X, E = load(src)
    cnt, own, hexes = {}, {}, []
    for nm, c in E.values():
        if nm in WALLS:
            if len(c) == 8:
                hexes.append(c)
            for h in HF:
                q = tuple(c[i] for i in h); k = frozenset(q)
                cnt[k] = cnt.get(k, 0) + 1; own[k] = q
    faces = [own[k] for k, v in cnt.items() if v == 1]
    wn = {n for q in faces for n in q}
    wed = [c for nm, c in E.values() if nm in SOL]
    assert wed, 'source needs build_batch312 solid lofts'
    # wedges near the attached row
    bad = {k for k, c in enumerate(wed) if any(v in wn for v in c)}
    for _ in range(ring):
        bn = {v for k in bad for v in wed[k]}
        bad |= {k for k, c in enumerate(wed) if any(v in bn for v in c)}
    keep = [c for k, c in enumerate(wed) if k not in bad]
    # outward faces of the kept wedges (penta6: bottom 0,1,2; top 3,4,5)
    WF = [((0, 2, 1), 'tri3'), ((3, 4, 5), 'tri3'), ((0, 1, 4, 3), 'quad4'), ((1, 2, 5, 4), 'quad4'), ((2, 0, 3, 5), 'quad4')]
    allwf = Counter(frozenset(c[i] for i in f) for c in wed for f, _ in WF)   # boundary of the WHOLE solid loft
    prim = [(g, [c[i] for i in f]) for c in keep for f, g in WF if allwf[frozenset(c[i] for i in f)] == 1]
    # nodes inside a wall hex at rest
    P = {v: X[v] for v in {v for _, f in prim for v in f}}
    tets = np.array([[X[h[i]] for i in tt] for h in hexes for tt in T5]); cen = tets.mean(1)
    tr = cKDTree(cen); rmax = np.max(np.linalg.norm(tets - cen[:, None], axis=2))
    inside = set()
    for v, p in P.items():
        for j in tr.query_ball_point(p, rmax):
            a, b, c_, d = tets[j]
            lam = np.linalg.solve(np.column_stack([b - a, c_ - a, d - a]), p - a)
            if min(1 - lam.sum(), lam.min()) > 1e-6:
                inside.add(v); break
    prim = [(g, f) for g, f in prim if not set(f) & inside]
    pn = np.array([X[v] for v in {v for _, f in prim for v in f}])
    lumen = set()
    for nm in ('SlidingElastic1Primary', 'SlidingElastic1Secondary'):
        a = t.index(f'<Surface name="{nm}"'); b = t.index('</Surface>', a)
        lumen |= {frozenset(int(v) for v in r.split(',')) for r in re.findall(r'>([0-9,]+)</', t[a:b])}
    ptr = cKDTree(pn)
    sec = [q for q in faces if frozenset(q) not in lumen and ptr.query(np.mean([X[v] for v in q], 0))[0] < dmax]
    # (L318b, sr 3) drop wall faces with a primary node BEHIND them at rest within SR (projecting inside the face): the
    # contact would kick those nodes at the first step (sr 0.5 never reached them)
    keep, behind = [], 0
    for q in sec:
        P = [X[v] for v in q]; cq = np.mean(P, 0); rad = max(np.linalg.norm(p - cq) for p in P)
        nq = np.cross(P[2] - P[0], P[3] - P[1]); nq /= np.linalg.norm(nq)
        jj = ptr.query_ball_point(cq, sr + rad)
        if jj:
            d = (pn[jj] - cq) @ nq; r = np.linalg.norm(pn[jj] - cq - d[:, None] * nq, axis=1)
            if np.any((d < 0) & (d > -sr) & (r < rad)):
                behind += 1; continue
        keep.append(q)
    sec = keep
    # smallest rest gap primary node -> secondary face centre (rough check against sr)
    gap = min(ptr.query(np.array([np.mean([X[v] for v in q], 0) for q in sec]))[0])
    surf = ('\t\t<Surface name="AVW_Para_solid_walls_primary">\n' +
            ''.join(f'\t\t\t<{g} id="{k}">{",".join(map(str, f))}</{g}>\n' for k, (g, f) in enumerate(prim, 1)) +
            '\t\t</Surface>\n\t\t<Surface name="AVW_Para_solid_walls_secondary">\n' +
            ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q))}</quad4>\n' for k, q in enumerate(sec, 1)) +
            '\t\t</Surface>\n\t\t<SurfacePair name="AVW_Para_solid_walls">\n\t\t\t<primary>AVW_Para_solid_walls_primary</primary>\n'
            '\t\t\t<secondary>AVW_Para_solid_walls_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', 'name="AVW_Para_solid_walls" surface_pair="AVW_Para_solid_walls"')
         .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
         .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
    for s_ in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>'):
        assert s_ in c, s_
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE: contact AVW_Para_solid_walls (sliding-elastic PENALTY {pen:g} auto, two_pass 0, search_radius {sr:g}, '
           f'offset 0): primary = {len(prim)} outward faces of the solid lofts ({len(bad)} of {len(wed)} wedges on / within {ring} '
           f'of the attached row left out; {len(inside)} nodes inside a wall hex at rest, their faces left out), secondary = '
           f'{len(sec)} wall boundary faces within {dmax:g} mm (lumen faces excluded; {behind} with a loft node behind them within {sr:g} mm at rest left out); nearest secondary face centre to a primary node at rest {gap:.2f} mm')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch318.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--sr', 0.5), g('--pen', 0.5), g('--dmax', 10.0), int(g('--ring', 1)))
