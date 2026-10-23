"""Batch 308 (2026-10-17; the user: "I wonder if a tighter mesh on the PM arch would help the contacts").
NOT IN SOURCE (mesh refinement): every PM_Plane quad split 2 x 2 (new mid-edge nodes shared between neighbours, one centre node
each; ~1-2 mm elements, not the 0.36-0.68 mm strips that pm_coarsen.py removed). Mid-edge nodes on edges with both ends in
BC-PM-outer-arc join that set (clamped). Every <Surface> made of PM_Plane facets (e.g. AVW_Para_PM_secondary,
AVW_Para_rope_PM_secondary) is rebuilt from the sub-quads with the same winding as its parent facet. Linear-constraint ties to
PM nodes are left as they are (the original corner nodes still exist).
usage: py -3.10 build_batch308.py SRC_RUN NEW_RUN"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new):
    t, X, E = load(src)
    pm = {e: c for e, (nm, c) in E.items() if nm == 'PM_Plane'}
    nid = max(X) + 1; eid = max(E) + 1
    mids, rows = {}, []
    def node(p):
        nonlocal nid
        rows.append(f'\t\t\t<node id="{nid}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>\n'); nid += 1
        return nid - 1
    def mid(a, b):
        k = frozenset((a, b))
        if k not in mids:
            mids[k] = node((X[a] + X[b]) / 2)
        return mids[k]
    sub = {}
    for e, (a, b, c, d) in pm.items():
        ab, bc, cd, da = mid(a, b), mid(b, c), mid(c, d), mid(d, a)
        m = node((X[a] + X[b] + X[c] + X[d]) / 4)
        sub[e] = [(a, ab, m, da), (ab, b, bc, m), (m, bc, c, cd), (da, m, cd, d)]
    # elements: replace the PM_Plane block
    blk = re.search(r'(\t\t<Elements type="quad4" name="PM_Plane">\n)(.*?)(\t\t</Elements>\n)', t, re.S)
    new_rows = []
    for e in pm:
        for k, q in enumerate(sub[e]):
            new_rows.append(f'\t\t\t<elem id="{e if k == 0 else eid}">{",".join(map(str, q))}</elem>\n')
            if k:
                eid += 1
    # the element ids above must stay increasing within the block: keep originals first, the new ones after
    first = [r for r in new_rows if int(re.search(r'id="(\d+)"', r).group(1)) in pm]
    rest = [r for r in new_rows if int(re.search(r'id="(\d+)"', r).group(1)) not in pm]
    t = t.replace(blk.group(0), blk.group(1) + ''.join(first + rest) + blk.group(3), 1)
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="PM_refine_nodes">\n' + ''.join(rows) + '\t\t</Nodes>\n' + t[i:]
    # BC set
    ns = re.search(r'<NodeSet name="BC-PM-outer-arc">([^<]*)</NodeSet>', t)
    bcn = [int(v) for v in re.findall(r'\d+', ns.group(1))]; sb = set(bcn)
    add = [v for k, v in mids.items() if k <= sb]
    t = t.replace(ns.group(0), '<NodeSet name="BC-PM-outer-arc">' + ','.join(map(str, bcn + add)) + '</NodeSet>', 1)
    # surfaces made of PM facets
    key = {frozenset(c): e for e, c in pm.items()}
    nsurf = []
    for m in re.finditer(r'(\t\t<Surface name="([^"]*)">\n)(.*?)(\t\t</Surface>\n)', t, re.S):
        facets = [list(map(int, r.split(','))) for r in re.findall(r'<quad4 id="\d+">([0-9,]+)</quad4>', m.group(3))]
        if not facets or not all(frozenset(f) in key for f in facets):
            continue
        out = []
        for f in facets:
            e = key[frozenset(f)]; c = pm[e]
            same = any(f == c[k:] + c[:k] for k in range(4))
            for q in sub[e]:
                out.append(list(q) if same else list(q)[::-1])
        body = ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q))}</quad4>\n' for k, q in enumerate(out, 1))
        t = t.replace(m.group(0), m.group(1) + body + m.group(4), 1)
        nsurf.append(f'{m.group(2)} ({len(facets)} -> {len(out)})')
    msg = (f'NOT IN SOURCE (mesh refinement): PM_Plane {len(pm)} quads split 2 x 2 -> {4 * len(pm)} ({len(rows)} new nodes); '
           f'{len(add)} new mid-edge nodes added to BC-PM-outer-arc; surfaces rebuilt: {", ".join(nsurf) or "none"}')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch308.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
