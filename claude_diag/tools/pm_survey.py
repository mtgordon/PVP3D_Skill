"""The perineal membrane (PM_Plane, the source's display body; a fixed rigid body in FEBio) against what could attach to
it (2026-09-29, the user's request: make the PM a real structure, outer edge fixed as OPAL325_PM_mid in
Normal_Generic_copy.inp, the pieces at its inner arc attached to it instead of anchored).
Prints: the PM's boundary loop; every fully fixed BC set's distance to the PM surface and boundary; how many nodes of
each tissue domain lie near the PM; and, per connector family anchored near the PM, where the anchors sit on the loop.
usage: py -3.10 pm_survey.py [RUN]   (default L87_springs_newline_rhoi0)"""
import os
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from febmodel import Feb, dist_to_surface  # noqa: E402
from paths import RUNS_DIR  # noqa: E402


def pm_geometry(f, dom='PM_Plane'):
    et, d = f.elem_blocks[dom]
    ec = Counter()
    for c in d.values():
        for i in range(len(c)):
            ec[tuple(sorted((c[i], c[(i + 1) % len(c)])))] += 1
    adj = defaultdict(list)
    for (a, b), k in ec.items():
        if k == 1:
            adj[a].append(b)
            adj[b].append(a)
    start = min(adj)
    loop, prev = [start], None
    while True:
        nx = [n for n in adj[loop[-1]] if n != prev]
        if not nx or nx[0] == start:
            break
        prev = loop[-1]
        loop.append(nx[0])
    tris = []
    for c in d.values():
        tris.append(c[:3])
        if len(c) == 4:
            tris.append([c[0], c[2], c[3]])
    return loop, tris


def surf_dist(f, pts, tris):
    return np.array([r[0] for r in dist_to_surface(f, pts, tris)]) if len(pts) else np.array([])


def main(run='L87_springs_newline_rhoi0'):
    f = Feb(os.path.join(RUNS_DIR, run, run + '.feb'))
    loop, tris = pm_geometry(f)
    P = np.array([f.nodes[n] for n in loop])
    print(f'PM_Plane: {len(f.elem_blocks["PM_Plane"][1])} elements, boundary loop {len(loop)} nodes, '
          f'perimeter {np.linalg.norm(np.diff(np.vstack([P, P[:1]]), axis=0), axis=1).sum():.0f} mm')
    for bc in f.root.find('Boundary'):
        ns = bc.get('node_set')
        ids = f.nodesets.get(ns, [])
        if not ids:
            continue
        pts = np.array([f.nodes[n] for n in ids])
        dd = surf_dist(f, pts, tris)
        dl = np.array([np.linalg.norm(P - p, axis=1).min() for p in pts])
        print(f'  BC {ns:30s} {len(ids):3d} nodes ({bc.get("type")}): to the PM surface min / median / max '
              f'{dd.min():6.2f} / {np.median(dd):6.2f} / {dd.max():6.2f} mm; to the nearest boundary node median {np.median(dl):6.2f}')
    for dom in ('_PickedSet347', '_PickedSet64', '_PickedSet346', '_PickedSet66', 'LA_PCMPRM', 'LA_PCM', 'LA_ICM'):
        ids = f.domain_nodes(dom)
        dd = surf_dist(f, np.array([f.nodes[n] for n in ids]), tris)
        print(f'  {dom:14s} {len(ids):5d} nodes: within 1 / 3 / 5 mm of the PM surface: {(dd < 1).sum()} / '
              f'{(dd < 3).sum()} / {(dd < 5).sum()}; closest {dd.min():.2f} mm')
    # every spring set with an end within 3 mm of the PM
    print('  spring sets with an end within 3 mm of the PM surface:')
    for nm, prs in f.discsets.items():
        ends = sorted({n for p in prs for n in p})
        dd = surf_dist(f, np.array([f.nodes[n] for n in ends]), tris)
        k = int((dd < 3).sum())
        if k:
            near = [n for n, x in zip(ends, dd) if x < 3]
            pos = sorted(int(np.argmin(np.linalg.norm(P - f.nodes[n], axis=1))) for n in near)
            print(f'    {nm:28s} {len(prs):3d} springs: {k} of their {len(ends)} end nodes within 3 mm (closest '
                  f'{dd.min():.2f} mm); nearest loop positions {pos[0]}-{pos[-1]} of {len(loop)}')
    return f, loop, tris


if __name__ == '__main__':
    main(*sys.argv[1:])
