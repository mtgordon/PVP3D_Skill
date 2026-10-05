"""How much of the LA's displacement is its supports moving? Splits the LA displacement at converged states into the
tied edge (the LA nodes held to the ATLA / posterior-arcus truss chains by the LA_truss_ties linear constraints) and
the rest, and gives the chain strains. usage: py -3.10 la_support_split.py RUN [T ...]"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
run = sys.argv[1]
want = [float(v) for v in sys.argv[2:]] or [1.0]
root = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot()
mesh = root.find('Mesh')
la = set()
for e in mesh.findall('Elements'):
    if e.get('name', '').startswith('LA_'):
        la |= {int(v) for el in e for v in el.text.split(',')}
tied = {}
for c in root.find('Constraints'):
    if c.get('name') != 'LA_truss_ties':
        continue
    for lc in c.findall('linear_constraint'):
        a, b = [int(n.get('id')) for n in lc.findall('node')]
        tied[a] = b
tied_la = sorted(n for n in tied if n in la)
chains = {ds.get('name'): [tuple(int(v) for v in d.text.split(',')) for d in ds]
          for ds in mesh.findall('DiscreteSet') if 'springs' in ds.get('name', '') and ('ATLA' in ds.get('name') or 'Arcus' in ds.get('name'))}
x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
idx = {int(n): i for i, n in enumerate(x.node_ids)}
T = np.array([s[0] for s in x.states])
L = np.array([idx[n] for n in sorted(la)])
Tt = np.array([idx[n] for n in tied_la])
free = np.array([idx[n] for n in sorted(la - set(tied_la))])
print(f'{run}: {len(la)} LA nodes, {len(tied_la)} tied to the chains ({len(set(tied.values()))} chain nodes)')
for tw in want:
    k = int(np.argmin(np.abs(T - tw)))
    u = x.var(k, 'displacement')
    m = lambda I: np.linalg.norm(u[I], axis=1)
    q = lambda a: f'{np.median(a):.1f} / {np.percentile(a, 90):.1f} / {a.max():.1f}'
    print(f' t {T[k]:.3f}: |u| median / p90 / max  all LA {q(m(L))}; tied edge {q(m(Tt))}; the rest {q(m(free))}')
    for nm, prs in chains.items():
        s = []
        for a, b in prs:
            X0 = np.linalg.norm(x.X[idx[a]] - x.X[idx[b]])
            s.append(np.linalg.norm(x.X[idx[a]] + u[idx[a]] - x.X[idx[b]] - u[idx[b]]) / X0 - 1)
        print(f'   {nm}: strain min {min(s):+.3f} median {np.median(s):+.3f} max {max(s):+.3f} ({len(s)} springs)')
