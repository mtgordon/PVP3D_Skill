import sys, collections
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from febmodel import Feb
f = Feb(sys.argv[1])
# linear-constraint ties: LA node <-> truss node
ties = collections.defaultdict(list)
for lc in f.root.find('Constraints').iter('linear_constraint'):
    ns = [(int(n.get('id')), n.get('bc'), float(n.text)) for n in lc.findall('node')]
    if ns[0][1] == 'x':
        a, b = ns[0][0], ns[1][0]
        ties[b].append(a); ties[a].append(b)
elem_of = collections.defaultdict(set)
for name, (et, d) in f.elem_blocks.items():
    for conn in d.values():
        for n in conn: elem_of[n].add(name)
disc_of = collections.defaultdict(list)
for name, pairs in f.discsets.items():
    for a, b in pairs:
        disc_of[a].append(name); disc_of[b].append(name)
bcn = collections.defaultdict(list)
for bc in f.section('Boundary'):
    for n in f.nodesets.get(bc.get('node_set'), []): bcn[n].append(bc.get('name'))
for chain in ('Posterior_Arcus_Left_springs', 'Posterior_Arcus_Right_springs', 'ATLA_Left_springs', 'ATLA_Right_springs'):
    pairs = f.discsets[chain]
    order = [pairs[0][0]] + [b for _, b in pairs]
    L = [np.linalg.norm(f.nodes[a] - f.nodes[b]) for a, b in pairs]
    print(f'== {chain}: {len(pairs)} springs, L0 min/mean/max {min(L):.2f}/{np.mean(L):.2f}/{max(L):.2f} mm, total {sum(L):.1f} mm')
    for n in order:
        tl = ties.get(n, [])
        tiedom = collections.Counter(f.elem_owner.get(0, '') for _ in [])
        doms = collections.Counter(d for t in tl for d in elem_of[t])
        other = [d for d in disc_of[n] if d != chain]
        print(f'  {n}: LA-tied {len(tl):2d} {dict(doms)} | elems {sorted(elem_of[n]) or "-"} | other springs {collections.Counter(other) or "-"} | BC {bcn.get(n, "-")}')
# P-arcus fans
for fan in ('P-arcus-L_fan', 'P-arcus-R_fan'):
    fn = set(f.domain_nodes(fan))
    shared = collections.Counter()
    for n in fn:
        for d in elem_of[n]:
            if d != fan: shared[d] += 1
        for d in disc_of[n]: shared['disc:' + d] += 1
        for b in bcn.get(n, []): shared['BC:' + b] += 1
    print(f'== {fan}: {len(f.elem_blocks[fan][1])} elems, {len(fn)} nodes; shared -> {dict(shared)}')
