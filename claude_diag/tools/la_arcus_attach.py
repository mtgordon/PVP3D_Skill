"""How the LA and the P-arcus fans are held to the posterior-arcus chains (read-only).

usage: py -3.10 la_arcus_attach.py MODEL.feb
Per chain: LA nodes tied to each chain node (linear constraints) and their distance; fan nodes tied (if any);
connector springs; and the LA's mesh-boundary nodes near the arcus that are NOT tied, with their distance.
"""
import sys
from collections import Counter, defaultdict

import numpy as np

from febmodel import Feb

f = Feb(sys.argv[1])
ties = defaultdict(list)          # chain node -> [(constraint name, other node)]
for c in f.root.find('Constraints'):
    for lc in c.iter('linear_constraint'):
        ns = [(int(n.get('id')), n.get('bc')) for n in lc.findall('node')]
        if ns[0][1] == 'x':
            ties[ns[1][0]].append((c.get('name'), ns[0][0]))
la_elems = [conn for name in ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri') for conn in f.elem_blocks[name][1].values()]
la_nodes = {n for c in la_elems for n in c}
edges = Counter()
for c in la_elems:
    for i in range(len(c)):
        edges[tuple(sorted((c[i], c[(i + 1) % len(c)])))] += 1
la_bnd = {n for e, k in edges.items() if k == 1 for n in e}
conn_springs = defaultdict(int)
for name, pairs in f.discsets.items():
    if name.startswith('Parcus_conn'):
        for a, b in pairs:
            conn_springs[a] += 1
tied_la = {o for lst in ties.values() for nm, o in lst if nm == 'LA_truss_ties'}
for side in ('Left', 'Right'):
    pairs = f.discsets[f'Posterior_Arcus_{side}_springs']
    order = [pairs[0][0]] + [b for _, b in pairs]
    P = np.array([f.nodes[n] for n in order])
    def dchain(x):
        return min(np.linalg.norm(P[i] + np.clip((x - P[i]) @ (P[i + 1] - P[i]) / ((P[i + 1] - P[i]) @ (P[i + 1] - P[i])), 0, 1)
                                  * (P[i + 1] - P[i]) - x) for i in range(len(P) - 1))
    print(f'== Posterior_Arcus_{side}: node | LA nodes tied (max dist mm) | fan nodes tied | connector springs')
    n_la = n_fan = 0
    for n in order:
        la = [o for nm, o in ties.get(n, []) if nm == 'LA_truss_ties']
        fan = [o for nm, o in ties.get(n, []) if nm == 'Parcus_fan_ties']
        n_la += len(la); n_fan += len(fan)
        dmax = max((np.linalg.norm(f.nodes[o] - f.nodes[n]) for o in la), default=0)
        print(f'  {n}: {len(la):2d} ({dmax:.1f}) | {len(fan)} | {conn_springs.get(n, 0)}')
    near = sorted((dchain(f.nodes[n]), n) for n in la_bnd)
    free_near = [(d, n) for d, n in near if d < 10 and n not in tied_la]
    tied_bnd = [(d, n) for d, n in near if d < 10 and n in tied_la]
    print(f'  total LA nodes tied {n_la}, fan nodes tied {n_fan}; LA boundary nodes within 10 mm of this chain: '
          f'{len(tied_bnd)} tied, {len(free_near)} free (free ones at ' + ', '.join(f'{d:.1f}' for d, _ in free_near[:25]) + ' mm)')
