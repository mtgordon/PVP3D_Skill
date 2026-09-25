"""Track connector attachment nodes and adjacent LA elements through a run's converged states.

usage: conn_watch.py run_name
"""
import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

from xplt import Xplt
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

run = sys.argv[1]
d = os.path.join(RUNS_DIR, run)
root = ET.parse(os.path.join(d, run + '.feb')).getroot()
mesh = root.find('Mesh')
pairs = []
for ds in mesh.findall('DiscreteSet'):
    if ds.get('name', '').startswith('LA_sphincter'):
        pairs += [(ds.get('name'), *map(int, e.text.split(','))) for e in ds.findall('delem')]
dom_names = [dd.get('name') for dd in root.find('MeshDomains')]

x = Xplt(os.path.join(d, run + '.xplt'))
idx = {int(n): i for i, n in enumerate(x.node_ids)}
for i, dd in enumerate(x.domains):
    dd['name'] = dom_names[i] if i < len(dom_names) else f'disc{i}'
la_doms = [i for i, dd in enumerate(x.domains) if dd['name'].startswith('LA_')]
# element adjacency for LA nodes
adj = defaultdict(list)
for di in la_doms:
    for k, conn in enumerate(x.domains[di]['conn']):
        for n in conn:
            adj[int(n)].append((di, k))

conv = [i for i, s in enumerate(x.states) if s[1] == 0]
sel = conv[-6:]
X0 = x.X
print('connector: LA node -> PeB node | elongation (mm) at last states | min J of LA elems around LA node')
for name, a, b in pairs:
    ia, ib = idx[a], idx[b]
    L0 = np.linalg.norm(X0[ia] - X0[ib])
    el, js = [], []
    for s in sel:
        u = x.var(s, 'displacement')
        L = np.linalg.norm((X0[ia] + u[ia]) - (X0[ib] + u[ib]))
        el.append(L - L0)
        J = x.var(s, 'relative volume')
        js.append(min(J[di + 1][k] for di, k in adj[ia]) if adj[ia] else np.nan)
    print(f'{name[-13:]:13s} {a:6d}->{b:6d} L0={L0:5.2f} dL=' + ' '.join('%6.3f' % v for v in el) +
          ' | minJ=' + ' '.join('%6.3f' % v for v in js))
print('times:', ' '.join('%.4f' % x.states[s][0] for s in sel))
