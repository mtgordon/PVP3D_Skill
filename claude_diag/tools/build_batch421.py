"""Batch 421 (2026-10-22; the user: "try adding contact between the vaginal wall (probably mostly the curves) and the PM arch. I
notice that it's penetrating that and maybe that interaction would push it inwards? ... It needs to be there long term anyways").
In L413b the curves pass through the PM arch (PM_Plane, quad4 shells) from t ~0.9 (rest z -13..-7, 1.0-1.8 mm), the AVW and
PVW from t 1.0 (1.3 / 3.4 mm by t 1.036; rest z -7..-3). NOT IN SOURCE (the source has no wall <-> PM contact):
sliding-elastic NAME, primary = the outward boundary faces of the canal (curves canal_wrap_hex_*, AVW _PickedSet347, PVW
_PickedSet64; one union boundary, lumen faces (SlidingElastic1*) left out) with their rest centroid in ZLO..ZHI, secondary =
the PM_Plane quads wound so their normal faces the canal (towards the canal's x-y centre at that level; shell contact is one-
sided), PM facets with a primary node behind them within SR at rest left out. Settings as walls_LA with penalty PEN, two_pass 0,
search_radius SR, offset 0.
usage: py -3.10 build_batch421.py SRC_RUN NEW_RUN [--name canal_PM] [--sr 1.0] [--pen 1] [--zlo -30] [--zhi 5]"""
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

CANAL = ('_PickedSet347', '_PickedSet64')


def main(src, new, name, sr, pen, zlo, zhi):
    t, X, E = load(src)
    X = {k: np.asarray(p, float) for k, p in dict(X).items()}
    cnt, own = Counter(), {}
    for nm, c in E.values():
        if (nm in CANAL or nm.startswith('canal_wrap_hex')) and len(c) == 8:
            for h in HF:
                q = tuple(c[i] for i in h); k = frozenset(q); cnt[k] += 1; own[k] = q
    lumen = set()
    for nm in ('SlidingElastic1Primary', 'SlidingElastic1Secondary'):
        a = t.index(f'<Surface name="{nm}"'); b = t.index('</Surface>', a)
        lumen |= {frozenset(int(v) for v in r.split(',')) for r in re.findall(r'>([0-9,]+)</', t[a:b])}
    prim = [own[k] for k, v in cnt.items() if v == 1 and k not in lumen
            and zlo <= np.mean([X[n][2] for n in own[k]]) <= zhi]
    pn_ids = sorted({v for q in prim for v in q}); pn = np.array([X[v] for v in pn_ids])
    # canal x-y centre by level, from all canal nodes
    cn = np.array([X[v] for nm, c in E.values() if (nm in CANAL or nm.startswith('canal_wrap_hex')) for v in c])
    pm = [list(c) for nm, c in E.values() if nm == 'PM_Plane']
    # consistent winding over the arch (propagate across shared edges), then the side facing the canal by majority vote
    byedge = {}
    for k, q in enumerate(pm):
        for i in range(4):
            byedge.setdefault(frozenset((q[i], q[(i + 1) % 4])), []).append(k)
    seen, stack = {0}, [0]
    while stack:
        k = stack.pop(); q = pm[k]
        for i in range(4):
            a, b = q[i], q[(i + 1) % 4]
            for j in byedge[frozenset((a, b))]:
                if j in seen:
                    continue
                r = pm[j]
                if any((r[m], r[(m + 1) % 4]) == (a, b) for m in range(4)):   # same direction: flip the neighbour
                    pm[j] = [r[0], r[3], r[2], r[1]]
                seen.add(j); stack.append(j)
    assert len(seen) == len(pm)
    votes = 0
    for q in pm:
        P = [X[v] for v in q]; cq = np.mean(P, 0); nq = np.cross(P[2] - P[0], P[3] - P[1])
        lv = cn[np.abs(cn[:, 2] - cq[2]) < 5]
        lv = lv if len(lv) else cn
        votes += np.sign((np.r_[lv[:, :2].mean(0), cq[2]] - cq) @ nq)
    flipped = 0
    if votes < 0:
        pm = [[q[0], q[3], q[2], q[1]] for q in pm]; flipped = len(pm)
    sec, behind = [], 0
    ptr = cKDTree(pn)
    for q in pm:
        P = [X[v] for v in q]; cq = np.mean(P, 0)
        nq = np.cross(P[2] - P[0], P[3] - P[1]); nq /= np.linalg.norm(nq)
        rad = max(np.linalg.norm(p - cq) for p in P)
        jj = ptr.query_ball_point(cq, sr + rad)
        if jj:
            d = (pn[jj] - cq) @ nq; r = np.linalg.norm(pn[jj] - cq - d[:, None] * nq, axis=1)
            if np.any((d < 0) & (d > -sr) & (r < rad)):
                behind += 1; continue
        sec.append(q)
    gap = min(ptr.query(np.array([np.mean([X[v] for v in q], 0) for q in sec]))[0])
    surf = (f'\t\t<Surface name="{name}_primary">\n' +
            ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, f))}</quad4>\n' for k, f in enumerate(prim, 1)) +
            f'\t\t</Surface>\n\t\t<Surface name="{name}_secondary">\n' +
            ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q))}</quad4>\n' for k, q in enumerate(sec, 1)) +
            f'\t\t</Surface>\n\t\t<SurfacePair name="{name}">\n\t\t\t<primary>{name}_primary</primary>\n'
            f'\t\t\t<secondary>{name}_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', f'name="{name}" surface_pair="{name}"')
    for key, val in (('penalty', f'{pen:g}'), ('two_pass', '0'), ('search_radius', f'{sr:g}'), ('offset', '0')):
        c = re.sub(r'<%s>[^<]*</%s>' % (key, key), f'<{key}>{val}</{key}>', c)
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE: contact {name} (sliding-elastic as walls_LA, penalty {pen:g}, two_pass 0, search_radius {sr:g}, '
           f'offset 0): primary = {len(prim)} outward canal boundary faces (curves + AVW + PVW, lumen faces out) with rest '
           f'centroid z {zlo:g}..{zhi:g}; secondary = {len(sec)} of {len(pm)} PM_Plane quads wound to face the canal ({flipped} '
           f'flipped after a consistent winding, vote {votes:+.0f}; {behind} with a primary node behind them within {sr:g} mm at rest left out); nearest PM face centre to a '
           f'primary node at rest {gap:.2f} mm')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch421.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], a[a.index('--name') + 1] if '--name' in a else 'canal_PM', g('--sr', 1.0), g('--pen', 1.0),
         g('--zlo', -30.0), g('--zhi', 5.0))
