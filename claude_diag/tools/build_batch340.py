"""Batch 340 (2026-10-19; the user: contact between the fan strap (build_batch337) and the curved edges). build_batch318 (solid
loft <-> wall contact) with the loft domains SOL, the target domains TGT and the contact NAME as arguments, and the attached row
found also by DISTANCE (a loft node within 0.05 mm of a target-surface node at rest, e.g. the spring-patch anchors the strap's
curve edge sits on), not only by shared nodes. NOT IN SOURCE: sliding-elastic NAME, primary = the outward faces of the solid
loft wedges away from the attached row (--ring), secondary = the target domains' boundary faces within DMAX mm, lumen faces
(SlidingElastic1*) and faces with a loft node behind them within SR at rest left out; penalty PEN, two_pass 0, search_radius
SR, offset 0, other settings as walls_LA.
usage: py -3.10 build_batch340.py SRC_RUN NEW_RUN --sol A,B --tgt C,D --name NAME [--sr 0.5] [--pen 0.5] [--dmax 10] [--ring 1]"""
import os
import re
import sys
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch297 import HF  # noqa: E402

T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]


def main(src, new, sr, pen, dmax, ring, SOL, WALLS, NAME):
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
    assert wed, ('no solid-loft wedges in', SOL)
    # wedges near the attached row
    wtree = cKDTree(np.array([X[v] for v in wn]))
    near = {v for c in wed for v in c if wtree.query(X[v])[0] < 0.05}
    bad = {k for k, c in enumerate(wed) if any(v in wn or v in near for v in c)}
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
    surf = (f'\t\t<Surface name="{NAME}_primary">\n' +
            ''.join(f'\t\t\t<{g} id="{k}">{",".join(map(str, f))}</{g}>\n' for k, (g, f) in enumerate(prim, 1)) +
            f'\t\t</Surface>\n\t\t<Surface name="{NAME}_secondary">\n' +
            ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q))}</quad4>\n' for k, q in enumerate(sec, 1)) +
            f'\t\t</Surface>\n\t\t<SurfacePair name="{NAME}">\n\t\t\t<primary>{NAME}_primary</primary>\n'
            f'\t\t\t<secondary>{NAME}_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', f'name="{NAME}" surface_pair="{NAME}"')
         .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
         .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
    for s_ in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>'):
        assert s_ in c, s_
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE: contact {NAME} (sliding-elastic PENALTY {pen:g} auto, two_pass 0, search_radius {sr:g}, '
           f'offset 0): primary = {len(prim)} outward faces of the solid lofts ({len(bad)} of {len(wed)} wedges on / within {ring} '
           f'of the attached row left out; {len(inside)} nodes inside a wall hex at rest, their faces left out), secondary = '
           f'{len(sec)} wall boundary faces within {dmax:g} mm (lumen faces excluded; {behind} with a loft node behind them within {sr:g} mm at rest left out); nearest secondary face centre to a primary node at rest {gap:.2f} mm')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch340.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    gs = lambda k: tuple(a[a.index(k) + 1].split(','))   # noqa: E731
    main(a[0], a[1], g('--sr', 0.5), g('--pen', 0.5), g('--dmax', 10.0), int(g('--ring', 1)), gs('--sol'), gs('--tgt'),
         a[a.index('--name') + 1])
