import sys
import numpy as np
from febmodel import Feb, dist_to_surface

f = Feb(sys.argv[1])
in_elems = f.nodes_in_elements()

# which nodes carry any BC
bc_fixed = {}
for bc in f.section('Boundary'):
    ns = bc.get('node_set')
    dofs = ''.join(k[0] for k in ('x_dof', 'y_dof', 'z_dof')
                   if bc.find(k) is not None and bc.find(k).text.strip() == '1')
    for n in f.nodesets.get(ns, []):
        bc_fixed.setdefault(n, set()).add((bc.get('name'), dofs))

# discrete connectivity
disc_deg = {}
for name, pairs in f.discsets.items():
    for a, b in pairs:
        disc_deg.setdefault(a, []).append(name)
        disc_deg.setdefault(b, []).append(name)

for bcname in ('BC-ATLA-L', 'BC-ATLA-R', 'BC-PosArcus_L', 'BC-PosArcus_R'):
    ids = f.nodesets[bcname]
    print(bcname, ids, [f.node_block.get(i) for i in ids])

print()
for cname in f.section('Contact'):
    sp = cname.get('surface_pair')
    if 'LA_ICM' not in sp:
        continue
    prim, sec = f.surfpairs[sp]
    md = float(cname.find('max_distance').text)
    tris = f.surface_tris(prim)
    snodes = f.surface_nodes(sec)
    d = dist_to_surface(f, [f.nodes[n] for n in snodes], tris)
    print(f'== {sp}: primary {prim} ({len(f.surfaces[prim])} facets), secondary {sec} '
          f'({len(f.surfaces[sec])} facets {set(t for t,_ in f.surfaces[sec])}, {len(snodes)} nodes), max_distance={md}')
    for n, (dist, ti) in zip(snodes, d):
        tag = []
        if n in bc_fixed: tag.append('FIXED:' + ','.join(sorted(x[0] + '/' + x[1] for x in bc_fixed[n])))
        if n not in in_elems: tag.append('no-elem')
        tag.append('springs=%d' % len(disc_deg.get(n, [])))
        eng = 'ENGAGE' if dist <= md else '      '
        print(f'  node {n:6d} {f.node_block[n]:32s} d={dist:7.3f} {eng} ' + ' '.join(tag))
