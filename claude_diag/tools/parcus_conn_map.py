"""Abaqus P-arcus-L/R-1..13 connectors (posterior-arcus truss -> VW-PeB) vs the FEBio model.

usage: py -3.10 parcus_conn_map.py MODEL.feb
For each connector: FEBio node of each end (by coordinate), initial length, angle to the chain tangent,
and whether the VW end is a node of the P-arcus fan. Then where the fan's loose (arcus-side) edge sits
relative to the chain.
"""
import re
import sys
from collections import Counter, defaultdict

import numpy as np

from abq_surf import lines
from febmodel import Feb
from inpmap import part_nodes, nearest_feb

f = Feb(sys.argv[1])
L = lines()
conns = []
for i, ln in enumerate(L):
    if ln.strip().lower().startswith('*element, type=conn3d2'):
        v = [t.strip() for t in L[i + 1].split(',')]
        m = re.search(r'elset=(P-arcus-[LR]-\d+)', L[i + 2])
        if m:
            conns.append((m.group(1), v[1], v[2]))
print(f'{len(conns)} P-arcus connectors in the .inp')

parts = {}
def pnodes(p):
    if p not in parts:
        parts[p] = part_nodes(p)
    return parts[p]

elem_of = defaultdict(set)
for name, (et, d) in f.elem_blocks.items():
    for conn in d.values():
        for n in conn:
            elem_of[n].add(name)

for side in 'LR':
    chain = f.discsets[f'Posterior_Arcus_{"Left" if side == "L" else "Right"}_springs']
    order = [chain[0][0]] + [b for _, b in chain]
    fan = f'P-arcus-{side}_fan'
    fan_nodes = set(f.domain_nodes(fan))
    print(f'\n== side {side}: connector | arcus node -> FEBio (d) | VW node -> FEBio (d, domains) | L0 | angle to chain')
    vw_feb = []
    for name, a, b in conns:
        if not name.startswith(f'P-arcus-{side}-'):
            continue
        pa, na = a.rsplit('.', 1)
        pb, nb = b.rsplit('.', 1)
        pa = pa.strip('"')[:-2]
        pb = pb.strip('"')[:-2]
        xa = pnodes(pa)[int(na)]
        xb = pnodes(pb)[int(nb)]
        fa, da = nearest_feb(f, xa)
        fb, db = nearest_feb(f, xb)
        vw_feb.append(fb)
        k = order.index(fa) if fa in order else -1
        if 0 < k < len(order) - 1:
            t = f.nodes[order[k + 1]] - f.nodes[order[k - 1]]
        elif k == 0:
            t = f.nodes[order[1]] - f.nodes[order[0]]
        else:
            t = None
        c = xb - xa
        ang = np.degrees(np.arccos(abs(c @ t) / np.linalg.norm(c) / np.linalg.norm(t))) if t is not None else float('nan')
        doms = ','.join(sorted(elem_of[fb]))
        infan = 'FAN' if fb in fan_nodes else '   '
        print(f'  {name:14s} {na:>3s} -> {fa} ({da:.1e}) | {nb:>4s} -> {fb} ({db:.1e}, {doms}) {infan} | {np.linalg.norm(c):5.2f} | {ang:5.1f} deg')
    # fan boundary: nodes on edges used by only one fan element
    ecount = Counter()
    for conn in f.elem_blocks[fan][1].values():
        for i in range(len(conn)):
            e = tuple(sorted((conn[i], conn[(i + 1) % len(conn)])))
            ecount[e] += 1
    bnodes = {n for e, c in ecount.items() if c == 1 for n in e}
    shared = {n for n in fan_nodes if len(elem_of[n]) > 1}
    loose = sorted(bnodes - shared)
    P = np.array([f.nodes[n] for n in order])
    def dist_to_chain(x):
        best = 1e9
        for i in range(len(P) - 1):
            a, b = P[i], P[i + 1]
            s = np.clip((x - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1)
            best = min(best, np.linalg.norm(a + s * (b - a) - x))
        return best
    dl = [dist_to_chain(f.nodes[n]) for n in loose]
    print(f'  fan {fan}: {len(fan_nodes)} nodes, {len(bnodes)} boundary, shared with other domains {len(shared)}'
          f' (VW connector ends among them: {len(set(vw_feb) & shared)}/{len(vw_feb)})')
    print(f'  loose boundary nodes: {len(loose)}; distance to the chain polyline [mm]: '
          + ' '.join(f'{d:.2f}' for d in sorted(dl)))
