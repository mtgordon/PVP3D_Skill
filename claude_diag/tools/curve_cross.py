"""Curve surface crossings (2026-10-22; the user on L399c: "how isn't the curve through itself here? ... something not correct about
how it's being looked at"). curve_selfpen.py counts only curve surface NODES that end inside another curve hex and gives their
distance to the nearest curve surface: for a 1-2 mm curve wall that depth tops out at about half the wall, and a node that has
gone right through and out the other side is not counted at all. This tool instead finds curve surface EDGES that pierce curve
surface triangles they do not touch (a true surface-surface crossing, whatever the depth), minus the crossings already present
at rest; per side the number of piercing edges, how far the piercing edges' nearer end has gone past the face (a lower bound of the
overlap), and the rest-position centroid of each crossing cluster.
usage: py -3.10 curve_cross.py RUN [--times 0.9,1.0,last]"""
import os
import sys
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
FI = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def surface(x):
    H = np.vstack([d['conn'][:, :8] for d in x.domains if (d.get('name') or '').startswith('canal_wrap_hex')])
    cnt = Counter(tuple(sorted(h[list(f)])) for h in H for f in FI)
    quads = [h[list(f)] for h in H for f in FI if cnt[tuple(sorted(h[list(f)]))] == 1]
    tris = np.array([q[[0, 1, 2]] for q in quads] + [q[[0, 2, 3]] for q in quads])
    edges = np.array(sorted({tuple(sorted((q[k], q[(k + 1) % 4]))) for q in quads for k in range(4)}))
    return tris, edges


def crossings(P, tris, edges):
    TP = P[tris]; cen = TP.mean(1); rad = np.max(np.linalg.norm(TP - cen[:, None], axis=2), 1)
    tr = cKDTree(cen); RM = rad.max()
    out = []
    for a, b in edges:
        A, B = P[a], P[b]; mid = (A + B) / 2; half = np.linalg.norm(B - A) / 2
        for j in tr.query_ball_point(mid, RM + half):
            t = tris[j]
            if a in t or b in t:
                continue
            p0, p1, p2 = TP[j]
            n = np.cross(p1 - p0, p2 - p0)
            da, db = (A - p0) @ n, (B - p0) @ n
            if da * db >= 0:
                continue
            X = A + (B - A) * (da / (da - db))
            v0, v1, v2 = p1 - p0, p2 - p0, X - p0
            d00, d01, d11, d20, d21 = v0 @ v0, v0 @ v1, v1 @ v1, v2 @ v0, v2 @ v1
            den = d00 * d11 - d01 * d01
            l1 = (d11 * d20 - d01 * d21) / den; l2 = (d00 * d21 - d01 * d20) / den
            if l1 >= 0 and l2 >= 0 and l1 + l2 <= 1:
                out.append((a, b, j, min(abs(da), abs(db)) / np.linalg.norm(n)))
    return out


def main(run, times):
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = converged_state_indices(x, os.path.join(RUNS, run, run + '.log'))
    tt = np.array([x.states[i][0] for i in idx])
    tris, edges = surface(x)
    X0 = x.X
    base = {(a, b, j) for a, b, j, d in crossings(X0, tris, edges)}
    print(f'{run}: {len(tris)} curve surface triangles, {len(edges)} edges; crossings at rest {len(base)}')
    for T in times:
        i = idx[-1] if T == 'last' else idx[int(np.argmin(abs(tt - float(T))))]
        P = X0 + x.var(i, 'displacement')
        cr = [c for c in crossings(P, tris, edges) if c[:3] not in base]
        side = Counter('x<0' if X0[a][0] < 0 else 'x>0' for a, b, j, d in cr)
        deep = [d for a, b, j, d in cr]
        pts = np.array([(X0[a] + X0[b]) / 2 for a, b, j, d in cr]) if cr else np.zeros((0, 3))
        clusters = []
        if len(pts):
            lab = -np.ones(len(pts), int); c = 0
            ptr = cKDTree(pts)
            for k in range(len(pts)):
                if lab[k] >= 0:
                    continue
                stack = [k]; lab[k] = c
                while stack:
                    m = stack.pop()
                    for n in ptr.query_ball_point(pts[m], 3.0):
                        if lab[n] < 0:
                            lab[n] = c; stack.append(n)
                c += 1
            dd = np.array(deep)
            clusters = sorted(((np.sum(lab == g), np.round(pts[lab == g].mean(0), 1), dd[lab == g].max()) for g in range(c)), key=lambda r: -r[0])
        print(f'  t {x.states[i][0]:.3f}: piercing edges {len(cr)} ({dict(side)}), of them > 0.5 mm deep {sum(d > 0.5 for d in deep)}, '
              f'> 1 mm {sum(d > 1 for d in deep)}; clusters (count, rest centroid, deepest edge end past the face): '
              + '; '.join(f'{n} at {tuple(cc)} {dm:.2f} mm' for n, cc, dm in clusters[:6]), flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    times = ['last']
    if '--times' in a:
        i = a.index('--times'); times = a[i + 1].split(','); a = a[:i] + a[i + 2:]
    main(a[0], times)
