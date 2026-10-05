"""Everything attached to a domain's nodes in a .feb: other domains sharing nodes, discrete springs (with where their
other end is), BCs, contact surfaces, surface loads and linear constraints.

usage: py -3.10 what_holds.py MODEL.feb [DOMAIN]     (default DOMAIN _PickedSet66, the perineal body)
"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from febmodel import Feb  # noqa: E402


def main():
    path = sys.argv[1]
    dom = sys.argv[2] if len(sys.argv) > 2 else '_PickedSet66'
    fe = Feb(path)
    root = fe.root
    P = set(fe.domain_nodes(dom))
    node_doms = defaultdict(set)
    for name, (et, d) in fe.elem_blocks.items():
        for conn in d.values():
            for n in conn:
                node_doms[n].add(name)
    # element blocks -> domain names (MeshDomains reference the Elements blocks by name)
    print(f'{os.path.basename(path)}: {dom}, {len(P)} nodes')
    print('-- other element blocks sharing nodes')
    c = Counter(b for n in P for b in node_doms[n] if b != dom)
    for b, k in c.most_common():
        print(f'   {k:5d} nodes shared with {b} ({fe.elem_blocks[b][0]})')
    # BCs
    bc_nodes = {}
    bnd = root.find('Boundary')
    fixed = defaultdict(set)
    if bnd is not None:
        for bc in bnd:
            ns = bc.get('node_set')
            ids = set(fe.nodesets.get(ns, []))
            dofs = ''.join(ax for ax in 'xyz' if (bc.findtext(f'{ax}_dof') or '0').strip() == '1')
            for n in ids:
                fixed[n].add(dofs)
            k = len(ids & P)
            if k:
                print(f'-- BC {bc.get("name")} ({dofs} fixed): {k} of its {len(ids)} nodes on {dom}')
    # discrete springs
    disc = root.find('Discrete')
    dmats = []
    if disc is not None:
        for dm in disc.findall('discrete_material'):
            dmats.append(dm.get('name'))
        print('-- discrete sets with an end on the domain')
        for d in disc.findall('discrete'):
            ds = d.get('discrete_set')
            pairs = fe.discsets.get(ds, [])
            on = [(a, b) for a, b in pairs if a in P or b in P]
            if not on:
                continue
            other = Counter()
            for a, b in on:
                o = b if a in P else a
                if o in P:
                    other['(both ends on the domain)'] += 1
                    continue
                where = sorted(node_doms[o]) or ['(no element)']
                fx = ''.join(sorted(fixed.get(o, set())))
                other[', '.join(where) + (f' [fixed {fx}]' if fx else '')] += 1
            mi = int(d.get('dmat')) - 1
            print(f'   {ds}: {len(on)} of {len(pairs)} springs (dmat {mi + 1} = {dmats[mi] if mi < len(dmats) else "?"}); '
                  f'other end: {dict(other)}')
    # contacts
    con = root.find('Contact')
    if con is not None:
        print('-- contacts')
        for ct in con:
            sp = ct.get('surface_pair')
            if sp not in fe.surfpairs:
                continue
            for role, sname in zip(('primary', 'secondary'), fe.surfpairs[sp]):
                facets = fe.surfaces.get(sname, [])
                k = sum(1 for _, conn in facets if all(n in P for n in conn))
                kn = len(set(fe.surface_nodes(sname)) & P) if facets else 0
                if k or kn:
                    print(f'   {ct.get("name")} ({ct.get("type")}) {role} {sname}: {k} of {len(facets)} facets, {kn} nodes on {dom}')
    # loads
    lds = root.find('Loads')
    if lds is not None:
        print('-- surface loads')
        for sl in lds.findall('surface_load'):
            sname = sl.get('surface')
            facets = fe.surfaces.get(sname, [])
            k = sum(1 for _, conn in facets if all(n in P for n in conn))
            if k:
                print(f'   {sl.get("name")} ({sl.get("type")}): {k} of {len(facets)} facets on {dom}')
    # linear constraints
    cons = root.find('Constraints')
    if cons is not None:
        print('-- linear constraints')
        for c in cons:
            k = 0
            for lc in c.iter('linear_constraint'):
                ids = [int(ch.get('node')) for ch in lc.iter() if ch.get('node')]
                if any(n in P for n in ids):
                    k += 1
            if k:
                print(f'   {c.get("name")}: {k} linear constraints touch {dom}')


if __name__ == '__main__':
    main()
