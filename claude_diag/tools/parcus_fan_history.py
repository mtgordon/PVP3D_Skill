"""How is the P-arcus fan's arcus-side edge held in each model version?

usage: py -3.10 parcus_fan_history.py A.feb B.feb ...
For P-arcus-L/R_fan: boundary nodes not shared with the PVW/PeB solids, and what holds them
(other element domains, BCs, discrete sets, contact surfaces, linear constraints).
"""
import sys
from collections import Counter, defaultdict

from febmodel import Feb

for path in sys.argv[1:]:
    f = Feb(path)
    elem_of = defaultdict(set)
    for name, (et, d) in f.elem_blocks.items():
        for conn in d.values():
            for n in conn:
                elem_of[n].add(name)
    bcs = defaultdict(set)
    bnd = f.section('Boundary')
    if bnd is not None:
        for bc in bnd:
            for n in f.nodesets.get(bc.get('node_set'), []):
                bcs[n].add(bc.get('name'))
    disc = defaultdict(set)
    for name, pairs in f.discsets.items():
        for a, b in pairs:
            disc[a].add(name if not name.startswith('stab_') else 'stab_*')
            disc[b].add(name if not name.startswith('stab_') else 'stab_*')
    used_surfs = set()
    con = f.section('Contact')
    if con is not None:
        for c in con:
            sp = c.get('surface_pair')
            if sp in f.surfpairs:
                used_surfs.update(f.surfpairs[sp])
    surf_of = defaultdict(set)
    for s in used_surfs:
        for n in f.surface_nodes(s):
            surf_of[n].add(s)
    lc = defaultdict(int)
    cons = f.section('Constraints')
    if cons is not None:
        for l in cons.iter('linear_constraint'):
            for nd in l.findall('node'):
                lc[int(nd.get('id'))] += 1
    print(f'== {path.split(chr(92))[-1].split("/")[-1]}')
    for fan in ('P-arcus-L_fan', 'P-arcus-R_fan'):
        if fan not in f.elem_blocks:
            print(f'  {fan}: absent')
            continue
        ec = Counter()
        for conn in f.elem_blocks[fan][1].values():
            for i in range(len(conn)):
                ec[tuple(sorted((conn[i], conn[(i + 1) % len(conn)])))] += 1
        bn = {n for e, c in ec.items() if c == 1 for n in e}
        solid = {n for n in bn if any(d.startswith('_PickedSet') for d in elem_of[n])}
        rest = bn - solid
        holders = Counter()
        free = 0
        for n in rest:
            h = [d for d in elem_of[n] if d != fan] + sorted(bcs[n]) + sorted(disc[n]) \
                + sorted(surf_of[n]) + (['linear constraint'] if lc[n] else [])
            for x in h:
                holders[x] += 1
            free += not h
        print(f'  {fan}: {len(bn)} boundary nodes, {len(solid)} on PVW/PeB; other {len(rest)}: '
              f'held by {dict(holders) or "-"}; completely free {free}')
