"""Strap vs curve crossings (2026-10-22; the user on L413b: "I still see a little poking through"). Like curve_cross.py, but between
the parts: (a) curve surface edges piercing the strap's exterior faces (PM_avw_bottom_*_sheet_solid / _sheet_wide, penta6 split
into triangles) and (b) strap exterior edges piercing the curve's exterior faces. Crossings present at rest are not counted
(the strap's tied edge sits on the curve rim). Depth = the nearer edge end's distance past the face plane (a lower bound).
Clusters by the edge's rest midpoint (3 mm linkage).
usage: py -3.10 strap_cross.py RUN [--times 0.9,1.0,last]"""
import os
import sys
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xplt_reader import Xplt, converged_state_indices  # noqa: E402
from curve_cross import surface, crossings, RUNS  # noqa: E402

PF = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]


def strap_surface(x):
    W = [d['conn'][:, :6] for d in x.domains if (d.get('name') or '').startswith('PM_avw_bottom_')]
    W = np.vstack(W)
    cnt = Counter(tuple(sorted(w[list(f)])) for w in W for f in PF)
    faces = [w[list(f)] for w in W for f in PF if cnt[tuple(sorted(w[list(f)]))] == 1]
    tris = [f[[0, 1, 2]] for f in faces] + [f[[0, 2, 3]] for f in faces if len(f) == 4]
    edges = sorted({tuple(sorted((f[k], f[(k + 1) % len(f)]))) for f in faces for k in range(len(f))})
    return np.array(tris), np.array(edges)


def clusters(X0, cr):
    if not cr:
        return []
    pts = np.array([(X0[a] + X0[b]) / 2 for a, b, j, d in cr]); dd = np.array([d for *_, d in cr])
    lab = -np.ones(len(pts), int); c = 0; tr = cKDTree(pts)
    for k in range(len(pts)):
        if lab[k] >= 0:
            continue
        st = [k]; lab[k] = c
        while st:
            m = st.pop()
            for n in tr.query_ball_point(pts[m], 3.0):
                if lab[n] < 0:
                    lab[n] = c; st.append(n)
        c += 1
    return sorted(((np.sum(lab == g), np.round(pts[lab == g].mean(0), 1), dd[lab == g].max()) for g in range(c)), key=lambda r: -r[0])


def main(run, times):
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = converged_state_indices(x, os.path.join(RUNS, run, run + '.log'))
    tt = np.array([x.states[i][0] for i in idx])
    ct, ce = surface(x)
    st, se = strap_surface(x)
    X0 = x.X
    b1 = {c[:3] for c in crossings(X0, st, ce)}
    b2 = {c[:3] for c in crossings(X0, ct, se)}
    print(f'{run}: strap {len(st)} exterior triangles / {len(se)} edges; curve {len(ct)} / {len(ce)}; at rest {len(b1)} + {len(b2)}')
    for T in times:
        i = idx[-1] if T == 'last' else idx[int(np.argmin(abs(tt - float(T))))]
        P = X0 + x.var(i, 'displacement')
        c1 = [c for c in crossings(P, st, ce) if c[:3] not in b1]
        c2 = [c for c in crossings(P, ct, se) if c[:3] not in b2]
        for lab, cr in (('curve edges through strap faces', c1), ('strap edges through curve faces', c2)):
            d = [c[3] for c in cr]
            print(f'  t {x.states[i][0]:.3f} {lab}: {len(cr)}, > 0.25 mm {sum(v > 0.25 for v in d)}, > 0.5 mm '
                  f'{sum(v > 0.5 for v in d)}, max {max(d) if d else 0:.2f} mm; clusters: '
                  + '; '.join(f'{n} at {tuple(cc)} {dm:.2f}' for n, cc, dm in clusters(X0, cr)[:4]), flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    times = ['last']
    if '--times' in a:
        i = a.index('--times'); times = a[i + 1].split(','); a = a[:i] + a[i + 2:]
    main(a[0], times)
