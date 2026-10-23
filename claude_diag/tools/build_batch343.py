"""Batch 343 (2026-10-19 night; the user's go, item 2 of NEXT_SESSION_PROMPT_2026-10-19b): the AVW-Para <-> PM contact on a model
whose lofts are ALREADY solid (the springs line's fibre lofts, build_batch337 --fam AVW-Para; build_batch312 / 320 make the solids
and the contact together). NOT IN SOURCE: sliding-elastic NAME, primary = every boundary face of the penta6 domains SOL (faces
used by one wedge only, wound outward from the wedge centre), secondary = the TARGET shell facets, wound as one to face the lofts
(majority vote of the nearest loft node at rest); facets with a loft node behind them at rest (within SR, inside FAC x the facet size sideways) left out
(build_batch320's L321 rule); settings as walls_LA with penalty PEN auto, two_pass 0, search_radius SR, offset 0 (L312b's).
usage: py -3.10 build_batch343.py SRC_RUN NEW_RUN [--sol AVW-Para-L_sheet_solid,AVW-Para-R_sheet_solid] [--target PM_Plane]
[--name AVW_Para_solid_PM] [--pen 1] [--sr 3] [--fac 0.8]"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

# penta6 faces (FEBio node order: bottom 0-1-2, top 3-4-5)
PF = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]


def main(src, new, sol, target, name, pen, sr, fac):
    t, X, E = load(src)
    faces, lnodes = [], set()
    for d in sol:
        els = [c for nm, c in E.values() if nm == d]
        assert els and all(len(c) == 6 for c in els), d
        cnt = Counter(frozenset(e[i] for i in f) for e in els for f in PF)
        for e in els:
            cen = np.mean([X[v] for v in e], 0)
            for f in PF:
                q = [e[i] for i in f]
                if cnt[frozenset(q)] > 1:
                    continue
                P = [X[v] for v in q]; fc = np.mean(P, 0)
                nq = np.cross(P[2] - P[0], P[3] - P[1]) if len(q) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
                if np.dot(nq, fc - cen) < 0:
                    q = q[::-1]
                faces.append(q)
            lnodes |= set(e)
    LN = np.array([X[v] for v in sorted(lnodes)])
    # the target shell is wound consistently: flip it as a whole by the majority vote of "faces the nearest loft node" (L312b
    # reversed all of PM_Plane; a per-facet vote turns 2 facets at the PM edge the other way)
    tq = [q for nm, q in E.values() if nm in target]
    vote = 0
    for q in tq:
        P = [X[v] for v in q]; cen = np.mean(P, 0)
        nq = np.cross(P[2] - P[0], P[3] - P[1]) if len(q) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
        vote += 1 if np.dot(nq, LN[int(np.argmin(np.linalg.norm(LN - cen, axis=1)))] - cen) < 0 else -1
    rev = vote > 0
    sec, dropped = [], 0
    flip = len(tq) if rev else 0
    for q in tq:
        q = q[::-1] if rev else q
        P = [X[v] for v in q]; cen = np.mean(P, 0)
        nq = np.cross(P[2] - P[0], P[3] - P[1]) if len(q) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
        nq = nq / np.linalg.norm(nq); rad = max(np.linalg.norm(p - cen) for p in P)
        dd = (LN - cen) @ nq; r = np.linalg.norm(LN - cen - dd[:, None] * nq, axis=1)
        if np.any((dd < 0.05) & (dd > -sr) & (r < fac * rad)):
            dropped += 1; continue
        sec.append(q)
    assert sec
    g = lambda q: 'quad4' if len(q) == 4 else 'tri3'   # noqa: E731
    surf = (f'\t\t<Surface name="{name}_primary">\n' +
            ''.join(f'\t\t\t<{g(q)} id="{k}">{",".join(map(str, q))}</{g(q)}>\n' for k, q in enumerate(faces, 1)) +
            f'\t\t</Surface>\n\t\t<Surface name="{name}_secondary">\n' +
            ''.join(f'\t\t\t<{g(q)} id="{k}">{",".join(map(str, q))}</{g(q)}>\n' for k, q in enumerate(sec, 1)) +
            f'\t\t</Surface>\n\t\t<SurfacePair name="{name}">\n\t\t\t<primary>{name}_primary</primary>\n'
            f'\t\t\t<secondary>{name}_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', f'name="{name}" surface_pair="{name}"')
         .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
         .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
    for s_ in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>'):
        assert s_ in c, s_
    t = t[:b] + '\t\t' + c + t[b:]
    tc = Counter(len(q) for q in faces)
    msg = (f'NOT IN SOURCE: contact {name} (tools/build_batch343.py): {len(faces)} boundary faces of the solid lofts {"+".join(sol)} '
           f'({tc[3]} tri, {tc[4]} quad, wound outward) vs {"+".join(target)} ({len(sec)} facets, {flip} rewound (as one) to face the lofts, '
           f'{dropped} left out: a loft node behind them at rest within {sr:g} mm / {fac:g} x the facet size) (sliding-elastic as '
           f'walls_LA, penalty {pen:g} auto, two_pass 0, search_radius {sr:g}, offset 0)')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch343.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    gs = lambda k, d: tuple(a[a.index(k) + 1].split(',')) if k in a else d   # noqa: E731
    main(a[0], a[1], gs('--sol', ('AVW-Para-L_sheet_solid', 'AVW-Para-R_sheet_solid')), gs('--target', ('PM_Plane',)),
         a[a.index('--name') + 1] if '--name' in a else 'AVW_Para_solid_PM', g('--pen', 1.0), g('--sr', 3.0), g('--fac', 0.8))
