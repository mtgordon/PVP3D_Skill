"""Batch 404 (2026-10-22; the user on L399c: "a node sticking through at 0.5888. What is the spacing of the mesh on the strap?
Do we need a finer mesh so that the contact will catch things better?" -> "build it"). The fan strap (PM_avw_bottom_left/right_
sheet_solid: 32 penta6 a side, in-plane edges median 2.9 mm, 1.0-5.8; the curve near its end ~1.0 mm) refined 4x in its plane:
every wedge (c0 c1 c2 | c3 c4 c5) is split at its edge midpoints into 4 wedges (one layer kept), shared mid-edge nodes, the
parent's winding kept; the parent's fibre direction (ElementData mat_axis) copied to its 4 children; the strap's element
blocks rewritten at the end of the mesh with new ids. Ties: a mid-edge node between two ATTACHED nodes (shared with another
domain, in a constraint, or in a NodeSet) follows its edge: added to the NodeSet when both ends are in it (the PM-end anchors),
else a linear constraint u_m = (u_a + u_b) / 2 (strap_refine_ties; as the model's other ties), so every tied edge moves as the
coarse edge did. Every surface facet that is a face of a strap wedge is replaced by its 4 (triangle) or 2 (quad) children,
winding kept. NOT IN SOURCE (the strap is ours).
usage: py -3.10 build_batch404.py SRC_RUN NEW_RUN"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402

STRAP = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
W3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def vol(X, c):
    return sum(np.dot(np.cross(X[c[b]] - X[c[a]], X[c[d]] - X[c[a]]), X[c[e]] - X[c[a]]) / 6 for a, b, d, e in W3)


def main(src, new):
    t, X, E = load(src)
    X = dict(X)
    dn = defaultdict(set)
    for nm, c in E.values():
        for v in c:
            dn[v].add(nm)
    c0, c1 = t.find('<Constraints>'), t.find('</Constraints>')
    cons = {int(v) for v in re.findall(r'<node id="(\d+)"', t[c0:c1])}
    nsets = {m.group(1): {int(v) for v in re.findall(r'\d+', m.group(2))}
             for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', t)}
    innset = set().union(*nsets.values())
    attached = lambda v: bool(dn[v] - set(STRAP)) or v in cons or v in innset   # noqa: E731
    nid = max(X) + 1; mid = {}; newnodes = []

    def m(a, b):
        nonlocal nid
        k = frozenset((a, b))
        if k not in mid:
            mid[k] = nid; X[nid] = (X[a] + X[b]) / 2; newnodes.append(nid); nid += 1
        return mid[k]
    eid = max(E) + 1
    blocks, children, facemap = {}, {}, {}
    for s in STRAP:
        par = [(k, c) for k, (nm, c) in E.items() if nm == s]
        out = []
        for k, c in par:
            a0, a1, a2, b0, b1, b2 = c
            m01, m12, m20 = m(a0, a1), m(a1, a2), m(a2, a0)
            M01, M12, M20 = m(b0, b1), m(b1, b2), m(b2, b0)
            ch = [[a0, m01, m20, b0, M01, M20], [m01, a1, m12, M01, b1, M12], [m20, m12, a2, M20, M12, b2],
                  [m01, m12, m20, M01, M12, M20]]
            assert all(vol(X, q) > 0 for q in ch), k
            for q in ch:
                out.append((eid, q)); eid += 1
            children[k] = ch
            # face children: bottom / top triangles and the three side quads, by the parent face's node set
            for tri, sub in (((a0, a1, a2), [(a0, m01, m20), (m01, a1, m12), (m20, m12, a2), (m01, m12, m20)]),
                             ((b0, b1, b2), [(b0, M01, M20), (M01, b1, M12), (M20, M12, b2), (M01, M12, M20)])):
                facemap[frozenset(tri)] = (tri, sub)
            for (p, q), (P, Q), mm, MM in (((a0, a1), (b0, b1), m01, M01), ((a1, a2), (b1, b2), m12, M12),
                                           ((a2, a0), (b2, b0), m20, M20)):
                facemap[frozenset((p, q, Q, P))] = ((p, q, Q, P), [(p, mm, MM, P), (mm, q, Q, MM)])
        blocks[s] = (par, out)
    # ties for mid-edge nodes between two attached nodes
    lc, nsadd = [], defaultdict(list)
    for k, v in mid.items():
        a, b = tuple(k)
        if not (attached(a) and attached(b)):
            continue
        both = [n for n, ns in nsets.items() if a in ns and b in ns]
        if both:
            for n in both:
                nsadd[n].append(v)
            continue
        for d in 'xyz':
            lc.append(f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{v}" bc="{d}">1</node>\n\t\t\t\t<node id="{a}" bc="{d}">-0.5</node>\n'
                      f'\t\t\t\t<node id="{b}" bc="{d}">-0.5</node>\n\t\t\t</linear_constraint>\n')
    # ---- write
    for s in STRAP:   # drop the old element blocks
        t = re.sub(rf'\t\t<Elements type="penta6" name="{s}">\n.*?\t\t</Elements>\n', '', t, flags=re.S)
    i = t.rindex('\t\t</Nodes>\n') + len('\t\t</Nodes>\n')
    t = t[:i] + ('\t\t<Nodes name="strap_refine_nodes">\n' + ''.join(
        f'\t\t\t<node id="{v}">{X[v][0]:.7g},{X[v][1]:.7g},{X[v][2]:.7g}</node>\n' for v in newnodes) + '\t\t</Nodes>\n') + t[i:]
    i = t.rindex('\t\t</Elements>\n') + len('\t\t</Elements>\n')
    t = t[:i] + ''.join(f'\t\t<Elements type="penta6" name="{s}">\n' + ''.join(
        f'\t\t\t<elem id="{k}">{",".join(map(str, q))}</elem>\n' for k, q in blocks[s][1]) + '\t\t</Elements>\n' for s in STRAP) + t[i:]
    for s in STRAP:   # fibre directions: each parent's line x4, lids renumbered
        a = t.index(f'<ElementData type="mat_axis" elem_set="{s}">'); b = t.index('</ElementData>', a)
        rows = re.findall(r'<elem lid="(\d+)">(.*?)</elem>', t[a:b])
        assert len(rows) == len(blocks[s][0])
        body = ''.join(f'\t\t\t<elem lid="{4 * j + r + 1}">{val}</elem>\n' for j, (lid, val) in enumerate(rows) for r in range(4))
        head = t[a:t.index('\n', a) + 1]
        t = t[:a] + head + body + '\t\t' + t[b:]
    for n, vs in nsadd.items():
        t = re.sub(rf'(<NodeSet name="{re.escape(n)}">)([^<]*)(</NodeSet>)', lambda mm: f'{mm.group(1)}{mm.group(2)},{",".join(map(str, vs))}{mm.group(3)}', t)
    stats = defaultdict(int)

    def surf_sub(ms):
        name, body = ms.group(1), ms.group(2)
        out = []
        for g, fid, nodes in re.findall(r'\t\t\t<(quad4|tri3) id="(\d+)">([0-9,]+)</\1>\n', body):
            q = [int(v) for v in nodes.split(',')]
            key = frozenset(q)
            if key in facemap:
                ref, sub = facemap[key]
                # keep the facet's own winding: same cyclic order as ref -> as is, else reversed
                n0 = len(q); r = list(ref)
                same = any(q == r[j:] + r[:j] for j in range(n0))
                for f in sub:
                    out.append(('tri3' if len(f) == 3 else 'quad4', list(f) if same else list(f)[::-1]))
                stats[name] += 1
            else:
                out.append((g, q))
        return (f'\t\t<Surface name="{name}">\n' + ''.join(f'\t\t\t<{g} id="{i}">{",".join(map(str, q))}</{g}>\n'
                                                       for i, (g, q) in enumerate(out, 1)) + '\t\t</Surface>\n')
    t = re.sub(r'\t\t<Surface name="([^"]+)">\n(.*?)\t\t</Surface>\n', surf_sub, t, flags=re.S)
    if lc:
        blk = ('\t\t<constraint name="strap_refine_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
               '\t\t\t<maxaug>10</maxaug>\n' + ''.join(lc) + '\t\t</constraint>\n')
        i = t.index('\t</Constraints>')
        t = t[:i] + blk + t[i:]
    nw = sum(len(b[0]) for b in blocks.values())
    msg = (f'NOT IN SOURCE (the user 10-22): the fan strap refined 4x in-plane: {nw} wedges -> {4 * nw} (one layer), '
           f'{len(newnodes)} mid-edge nodes; fibre directions copied to the children; {len(lc) // 3} mid-edge nodes on tied edges '
           f'tied to their edge (strap_refine_ties), ' + ', '.join(f'{len(v)} added to NodeSet {n}' for n, v in nsadd.items()) +
           '; surface facets refined: ' + ', '.join(f'{n} {k}' for n, k in stats.items()))
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch404.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
