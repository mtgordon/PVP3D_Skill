"""Batch 302 (2026-10-17 night; the user away, loft-contacts plan step 1).
L296_chainsfan_paraPM_MAX (AVW-Para <-> PM arch contact, penalty 0.5) holds the AVW-Para sheets out of the PM but crawls from
t 0.3: the sheet's chain nodes slide over the PM arch's FREE inner edge (PM nodes 24191 / 24192 / 24198 are free boundary
nodes) and flip on and off its last facet (skill item 13, the contact-only seam stall; the user's ledges fixed that one).
NOT IN SOURCE (a convergence aid, like the ledges): a PM ledge = one row of quad4 shells extruded LEDGE mm in the PM's plane
outward from every PM boundary edge that is not on the clamped outer arc (BC-PM-outer-arc), sharing the edge nodes; PM_Yeoh,
thickness 2 (as PM_Plane); its facets are added to the AVW_Para_PM secondary surface (reversed like the PM's). Unloaded,
attached to nothing else.
usage: py -3.10 build_batch302.py SRC_RUN NEW_RUN [--ledge 3]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, L):
    t, X, E = load(src)
    pm = [c for nm, c in E.values() if nm == 'PM_Plane']
    ns = re.search(r'<NodeSet name="BC-PM-outer-arc">(.*?)</NodeSet>', t, re.S).group(1)
    fixed = {int(v) for v in re.findall(r'\d+', ns)}
    cnt = Counter(frozenset((q[i], q[(i + 1) % 4])) for q in pm for i in range(4))
    edges = []   # (a, b, owner quad) with a -> b in the owner's winding
    for q in pm:
        for i in range(4):
            a, b = q[i], q[(i + 1) % 4]
            if cnt[frozenset((a, b))] == 1 and not (a in fixed and b in fixed):
                edges.append((a, b, q))
    out = defaultdict(list)
    for a, b, q in edges:
        P = np.array([X[v] for v in q]); n = np.cross(P[2] - P[0], P[3] - P[1]); n /= np.linalg.norm(n)
        o = np.cross(X[b] - X[a], n); o /= np.linalg.norm(o)
        if np.dot(o, (X[a] + X[b]) / 2 - P.mean(0)) < 0:
            o = -o
        out[a].append(o); out[b].append(o)
    nid = max(X) + 1; eid = max(E) + 1
    new_of, rows_n = {}, []
    for v in sorted(out):
        d = np.mean(out[v], 0); d /= np.linalg.norm(d)
        p = X[v] + L * d; new_of[v] = nid
        rows_n.append(f'\t\t\t<node id="{nid}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>\n'); nid += 1
    quads = [(b, a, new_of[a], new_of[b]) for a, b, q in edges]   # owner traverses a->b, the neighbour b->a
    rows_e = [f'\t\t\t<elem id="{eid + k}">{",".join(map(str, c))}</elem>\n' for k, c in enumerate(quads)]
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="PM_ledge_nodes">\n' + ''.join(rows_n) + '\t\t</Nodes>\n' + t[i:]
    i = t.rindex('</Elements>') + len('</Elements>\n')
    t = t[:i] + '\t\t<Elements type="quad4" name="PM_ledge">\n' + ''.join(rows_e) + '\t\t</Elements>\n' + t[i:]
    dom = re.search(r'\t\t<ShellDomain name="PM_Plane".*?</ShellDomain>\n', t, re.S).group(0)
    t = t.replace(dom, dom + dom.replace('name="PM_Plane"', 'name="PM_ledge"'), 1)
    a = t.index('<Surface name="AVW_Para_PM_secondary">'); b = t.index('\t\t</Surface>', a)
    k0 = len(re.findall(r'<quad4 id=', t[a:b]))
    add = ''.join(f'\t\t\t<quad4 id="{k0 + k + 1}">{",".join(map(str, c[::-1]))}</quad4>\n' for k, c in enumerate(quads))
    t = t[:b] + add + t[b:]
    msg = (f'NOT IN SOURCE (convergence aid): PM ledge, {len(quads)} quad4 shells (PM_ledge, PM_Yeoh, thickness 2) extruded '
           f'{L:g} mm in the PM plane from the {len(edges)} PM boundary edges off the clamped outer arc ({len(rows_n)} new nodes '
           f'{min(new_of.values())}-{nid - 1}); added to the AVW_Para_PM secondary surface (reversed like the PM facets)')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch302.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--ledge') + 1]) if '--ledge' in a else 3.0)
