"""Audit which boundary conditions actually act on the model, across file versions.

For every <bc> in every .feb given, counts how many of its node_set's nodes are
"active" (used by an <Elements> block or a <DiscreteSet>) versus declared. A BC
whose nodes are all orphans pins nothing, and FEBio gives no warning (the orphans
just vanish in "N isolated vertices removed").

The classic cause is mesh surgery that renumbers nodes: a remesher (e.g. mmgs)
gives a fan/patch's boundary nodes new IDs at the same coordinates, the elements
move to the new IDs, and the BC node set keeps pointing at the old ones. For each
BC with orphan nodes, the audit lists live nodes at the *identical* coordinate
(within --tol) and which domain uses them: those are what the BC was meant to pin.

usage: py -3 bc_audit.py old.feb [newer.feb ...] [--tol 1e-5]
Needs Python 3 + numpy; imports feb_model.py from this directory.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feb_model import Feb  # noqa: E402


def audit(path, tol):
    f = Feb(path)
    active = set()
    node_dom = defaultdict(set)
    for name, (et, d) in f.elem_blocks.items():
        for conn in d.values():
            active.update(conn)
            for n in conn:
                node_dom[n].add(name)
    for pairs in f.discsets.values():
        for a, b in pairs:
            active.update((a, b))
    live = np.array(sorted(active))
    L = np.array([f.nodes[n] for n in live])
    rows = {}
    for bc in f.section('Boundary'):
        ns = bc.get('node_set')
        ids = f.nodesets.get(ns, [])
        on = [n for n in ids if n in active]
        hint = defaultdict(int)
        for n in ids:
            if n in active or n not in f.nodes:
                continue
            d = np.linalg.norm(L - f.nodes[n], axis=1)
            for k in np.where(d < tol)[0]:
                if int(live[k]) in ids:          # already pinned by this BC
                    continue
                for dom in node_dom[int(live[k])] or {'(discrete only)'}:
                    hint[dom] += 1
        rows[bc.get('name')] = (len(on), len(ids), dict(hint))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('febs', nargs='+')
    ap.add_argument('--tol', type=float, default=1e-5)
    a = ap.parse_args()
    res = [(p, audit(p, a.tol)) for p in a.febs]
    names = []
    for _, rows in res:
        names += [n for n in rows if n not in names]
    w = max(12, *(len(os.path.basename(p)[:28]) for p, _ in res))
    print('%-32s' % 'bc  (active/declared nodes)' + ''.join(('%' + str(w + 2) + 's') % os.path.basename(p)[:28] for p, _ in res))
    for n in names:
        print('%-32s' % n[:32] + ''.join(('%' + str(w + 2) + 's') % (
            '%d/%d' % rows[n][:2] if n in rows else '-') for _, rows in res))
    print()
    for p, rows in res:
        for n, (on, tot, hint) in rows.items():
            if on < tot and hint:
                print(f'{os.path.basename(p)}: {n}: {tot - on} orphan node(s); live nodes at the same '
                      f'coordinates are used by {hint} -> probably what this BC should pin')


if __name__ == '__main__':
    main()
