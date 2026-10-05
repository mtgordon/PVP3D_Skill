"""Deformed distance between the PVW_LA contact surfaces (the contact removed in the run): how close does the LA get to the
posterior vaginal wall? Also which surface is which."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from febmodel import Feb, pt_tri_dist
from xplt import Xplt
run = sys.argv[1]
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
f = Feb(os.path.join(R, run, run + '.feb'))
x = Xplt(os.path.join(R, run, run + '.xplt'))
idx = {int(n): i for i, n in enumerate(x.node_ids)}
conv = [i for i, s in enumerate(x.states) if s[1] == 0]
dom_of = {}
for nm, (t, d) in f.elem_blocks.items():
    for c in d.values():
        for n in c:
            dom_of.setdefault(n, set()).add(nm)
def tris(name):
    out = []
    for fac in f.surface_tris(name):
        out.append(fac)
    return out
P_tris, S_tris = f.surface_tris('PVW_LA_primary'), f.surface_tris('PVW_LA_secondary')
for nm, T in (('primary', P_tris), ('secondary', S_tris)):
    ns = {n for t in T for n in t}
    doms = {}
    for n in ns:
        for d in dom_of.get(n, ()):
            doms[d] = doms.get(d, 0) + 1
    print(nm, len(T), 'tris,', len(ns), 'nodes; domains', sorted(doms.items(), key=lambda kv: -kv[1])[:3])
for k in (conv[0], conv[len(conv)//2], conv[int(len(conv)*0.9)], conv[-1]):
    t = x.states[k][0]
    u = x.var(k, 'displacement')
    X = lambda n: x.X[idx[n]] + u[idx[n]]
    pn = sorted({n for tr in P_tris for n in tr})
    pts = np.array([X(n) for n in pn])
    A = np.array([X(tr[0]) for tr in S_tris]); B = np.array([X(tr[1]) for tr in S_tris]); C = np.array([X(tr[2]) for tr in S_tris])
    cen = (A + B + C) / 3
    d = []
    for p in pts:
        ii = np.argsort(np.linalg.norm(cen - p, axis=1))[:40]
        d.append(min(pt_tri_dist(p, A[i], B[i], C[i]) for i in ii))
    d = np.array(d)
    print(f't {t:.3f}: primary-node to secondary-surface distance min {d.min():.2f} mm, 1st pct {np.percentile(d, 1):.2f}; nodes < 2.0 mm (the offset): {(d < 2.0).sum()}, < 1.0: {(d < 1.0).sum()}, < 0.25: {(d < 0.25).sum()} of {len(d)}')
