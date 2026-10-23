"""Curve self-penetration (2026-10-14; the user: L283_junction_fine2z18_MAX "the sides still penetrate themselves around 0.65 on
one side and 0.69 on the other"). For each run and time: the curves' surface nodes (canal_wrap_hex_* boundary) that lie inside
another curve hex (hex8 split into 5 tets; the node's own hexes skipped), minus those already inside at rest; per side (x < 0 /
x > 0): the count and the deepest one (distance to the nearest curve boundary triangle not touching the node).
usage: py -3.10 curve_selfpen.py RUN [RUN ...] [--times 0.5,0.6,...]   (default: every 0.05 up to the last converged state)"""
import os
import sys
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt, converged_state_indices  # noqa: E402
from feb_model import pt_tri_dist  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
SHORT = {'canal_wrap_hex_avw': 'curveA', 'canal_wrap_hex_avw_cervix': 'curveAcx', 'canal_wrap_hex_pvw': 'curveP',
         'canal_wrap_hex_pvw_fine': 'curveP_fine'}
FI = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]


class Curves:
    def __init__(self, x):
        H, D = [], []
        for d in x.domains:      # every curve domain (2026-10-15: also _distal, _around, ...; the cervix part excluded)
            n = d.get('name') or ''
            if n.startswith('canal_wrap_hex') and n != 'canal_wrap_hex_avw_cervix' or n in SHORT:
                H.append(d['conn'][:, :8]); D += [SHORT.get(n, n.replace('canal_wrap_hex_', 'curve_'))] * len(d['conn'])
        self.H = H = np.vstack(H); self.D = np.array(D)
        cnt = Counter(tuple(sorted(h[list(f)])) for h in H for f in FI)
        bf = []
        for h in H:
            for f in FI:
                if cnt[tuple(sorted(h[list(f)]))] == 1:
                    q = h[list(f)]; bf += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        self.bf = np.array(bf); self.bnodes = np.unique(self.bf)
        self.tets = np.array([[h[i] for i in t] for h in H for t in T5]); self.tel = np.repeat(np.arange(len(H)), 5)
        self.elnodes = [set(h) for h in H]

    def inside(self, X):
        P = X[self.tets]; c = P.mean(1); r = np.max(np.linalg.norm(P - c[:, None], axis=2), 1)
        tree = cKDTree(c); R = r.max(); out = {}
        for n in self.bnodes:
            p = X[n]
            for t in tree.query_ball_point(p, R):
                if n in self.elnodes[self.tel[t]]:
                    continue
                a, b, cc, dd = P[t]
                try:
                    lam = np.linalg.solve(np.column_stack([b - a, cc - a, dd - a]), p - a)
                except np.linalg.LinAlgError:
                    continue
                if min(1 - lam.sum(), lam.min()) > 1e-3:
                    out[n] = self.tel[t]; break
        return out


def report(run, times=None):
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = converged_state_indices(x, os.path.join(RUNS, run, run + '.log'))
    tt = np.array([x.states[i][0] for i in idx])
    if times is None:
        times = [t for t in np.arange(0.05, 1.0001, 0.05) if t <= tt.max() + 1e-6] + [float(tt.max())]
    cv = Curves(x); base = set(cv.inside(x.X))
    print(f'{run}: last converged t {tt.max():.3f}; curve surface nodes inside another curve hex at rest: {len(base)}')
    seen = set()
    for T in times:
        i = idx[int(np.argmin(abs(tt - T)))]
        if i in seen:
            continue
        seen.add(i)
        X = x.X + x.var(i, 'displacement')
        hits = [n for n in cv.inside(X) if n not in base]
        tr = cKDTree(X[cv.bf].mean(1))
        side = {'x<0': [], 'x>0': []}
        for n in hits:
            best = 1e9
            for k in tr.query_ball_point(X[n], 6.0):
                if n not in cv.bf[k]:
                    best = min(best, pt_tri_dist(X[n], *X[cv.bf[k]]))
            side['x<0' if X[n, 0] < 0 else 'x>0'].append(best)
        s = '  '.join(f'{k}: {len(v)} nodes, deepest {max(v):.2f} mm' if v else f'{k}: 0' for k, v in side.items())
        print(f'  t {x.states[i][0]:.3f}  {s}', flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]; times = None
    if '--times' in a:
        i = a.index('--times'); times = [float(v) for v in a[i + 1].split(',')]; a = a[:i] + a[i + 2:]
    for r in a:
        report(r, times)
