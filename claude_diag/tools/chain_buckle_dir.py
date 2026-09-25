"""Which way does a compressed posterior-arcus chain node drift, relative to the (missing) P-arcus
connector direction and the LA surface normal?

usage: py -3.10 chain_buckle_dir.py RUN [L|R]
At the last converged state: each interior chain node's zig-zag drift w (its displacement minus the mean
of its two neighbours', across the chain), split into components along c_perp (the connector direction to the PVW, made
perpendicular to the chain) and along n_LA (mean normal of the LA elements tied to that node).
"""
import os
import sys
from collections import defaultdict

import numpy as np

from febmodel import Feb
from xplt import Xplt
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

RUNS = RUNS_DIR
# Abaqus P-arcus connector ends, FEBio ids (parcus_conn_map.py): chain node 23700+k (L) / 23719+k (R)
PVW = {'L': [1847, 1913, 1912, 1911, 1910, 1909, 1908, 1907, 1906, 1905, 1904, 1902, 1900],
       'R': [1865, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1980, 1982]}

run = sys.argv[1]
side = sys.argv[2] if len(sys.argv) > 2 else 'R'
d = os.path.join(RUNS, run)
f = Feb(os.path.join(d, run + '.feb'))
x = Xplt(os.path.join(d, run + '.xplt'))
idx = {int(n): i for i, n in enumerate(x.node_ids)}
s = [i for i, st in enumerate(x.states) if st[1] == 0][-1]
u = x.var(s, 'displacement')
print(f'{run}: last converged t = {x.states[s][0]:.4f}, side {side}')

pairs = f.discsets[f'Posterior_Arcus_{"Left" if side == "L" else "Right"}_springs']
order = [pairs[0][0]] + [b for _, b in pairs]
base = 23700 if side == 'L' else 23719
conn = {base + k + 1: PVW[side][k] for k in range(13)}

ties = defaultdict(list)
for lc in f.root.find('Constraints').iter('linear_constraint'):
    ns = [(int(n.get('id')), n.get('bc')) for n in lc.findall('node')]
    if ns[0][1] == 'x':
        ties[ns[1][0]].append(ns[0][0])
la_elems = defaultdict(list)
for name in ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri'):
    for conn_ in f.elem_blocks[name][1].values():
        for n in conn_:
            la_elems[n].append(conn_)

def pos(n):
    return f.nodes[n] + u[idx[n]]

def la_normal(nodes):
    acc = np.zeros(3)
    for n in nodes:
        for c in la_elems[n]:
            p = [pos(m) for m in c]
            nn = np.cross(p[1] - p[0], p[2] - p[0])
            if len(c) == 4:
                nn = nn + np.cross(p[2] - p[0], p[3] - p[0])
            acc += nn / np.linalg.norm(nn)
    return acc / np.linalg.norm(acc) if np.linalg.norm(acc) else acc

print(' node  drift mm | along c_perp | out of (t, c) plane | along n_LA | in LA plane | angle(c_perp, n_LA)')
for k in range(1, len(order) - 1):
    a, n, b = order[k - 1], order[k], order[k + 1]
    t = pos(b) - pos(a)
    t /= np.linalg.norm(t)
    # zig-zag drift: this node's displacement relative to its neighbours' mean, across the chain
    w = u[idx[n]] - 0.5 * (u[idx[a]] + u[idx[b]])
    w = w - (w @ t) * t
    nl = la_normal(ties[n])
    nl = nl - (nl @ t) * t
    nl /= np.linalg.norm(nl)
    if n in conn:
        c = pos(conn[n]) - pos(n)
        c = c - (c @ t) * t
        c /= np.linalg.norm(c)
        wc, ang = w @ c, np.degrees(np.arccos(abs(c @ nl)))
        wb = w @ np.cross(t, c)   # out of the plane of the chain and the connector (the fan's plane)
    else:
        wc, ang, wb = float('nan'), float('nan'), float('nan')
    wn = w @ nl
    rest = np.sqrt(max(w @ w - wn ** 2, 0))
    print(f'{n:6d} {np.linalg.norm(w):6.2f} | {wc:+7.2f}      | {wb:+7.2f}            | {wn:+7.2f}    | {rest:5.2f} | {ang:5.1f}'
          + ('   <- 23723' if n == 23723 else ''))
