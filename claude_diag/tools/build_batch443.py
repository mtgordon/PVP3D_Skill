"""Batch 443 (2026-10-23; the user on L440 / L441: "I don't think you attached the strap to the finer distal mesh correctly").
After the distal curve split around the U (build_batch286), every curve rim edge the strap is tied along (directly, shared
nodes) got a new curve node at its middle that is FREE, while the strap's own mid-edge node there (build_batch404, tied to the
straight line between the two rim nodes by strap_refine_ties) sits 0.0-0.22 mm away: strap and curve were joined only at the
coarse rim nodes. NOT IN SOURCE (our strap): each such strap node (a strap_refine_ties dependent whose hosts are both strap-curve
shared nodes) is replaced in the strap elements and every surface by the nearest curve node within TOL mm (the curve's rim
mid-node), and its strap_refine_ties rows are removed; the strap now shares the curve's fine rim nodes. Strap wedges checked
positive afterwards; the old strap node is left unused.
usage: py -3.10 build_batch443.py SRC_RUN NEW_RUN [--tol 0.3]"""
import os
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

W3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def main(src, new, tol):
    t, X, E = load(src)
    X = {k: np.asarray(v, float) for k, v in dict(X).items()}
    dn = defaultdict(set)
    for k, (n, c) in E.items():
        for v in c:
            dn[v].add(n)
    strap = {v for v, s in dn.items() if any(q.startswith('PM_avw_bottom') for q in s)}
    curve = {v for v, s in dn.items() if any(q.startswith('canal_wrap') for q in s)}
    shared = strap & curve
    blk = re.search(r'\t\t<constraint name="strap_refine_ties".*?</constraint>\n', t, re.S).group(0)
    rows = re.findall(r'\t\t\t<linear_constraint>.*?</linear_constraint>\n', blk, re.S)
    hosts = defaultdict(set)
    for r in rows:
        nn = [int(a) for a in re.findall(r'<node id="(\d+)" bc="\w">', r)]
        hosts[nn[0]] |= set(nn[1:])
    cfree = [v for v in curve if v not in strap]
    ct = cKDTree(np.array([X[v] for v in cfree]))
    merge = {}
    for s, h in hosts.items():
        if s in strap and h and h <= shared:
            d, j = ct.query(X[s])
            if d <= tol:
                merge[s] = (cfree[j], d)
    assert merge, 'nothing to merge'
    # strap elements
    for k, (n, c) in E.items():
        if n.startswith('PM_avw_bottom') and set(c) & set(merge):
            old = f'<elem id="{k}">' + ','.join(map(str, c)) + '</elem>'
            assert t.count(old) == 1, k
            t = t.replace(old, f'<elem id="{k}">' + ','.join(str(merge[v][0]) if v in merge else str(v) for v in c) + '</elem>')
    # surfaces
    def surf(m):
        return re.sub(r'(<(?:quad4|tri3) id="\d+">)([\d,]+)(</)',
                      lambda q: q.group(1) + ','.join(str(merge[int(v)][0]) if int(v) in merge else v for v in q.group(2).split(',')) + q.group(3),
                      m.group(0))
    t = re.sub(r'<Surface name="[^"]+">.*?</Surface>', surf, t, flags=re.S)
    # drop the merged nodes' tie rows
    keep = [r for r in rows if int(re.search(r'<node id="(\d+)"', r).group(1)) not in merge]
    head = blk[:blk.index('\t\t\t<linear_constraint>')]
    # an empty constraint block fails FEBio's read ("needs to have property linear_constraint"): drop it then
    t = t.replace(blk, head + ''.join(keep) + '\t\t</constraint>\n' if keep else '')
    # check the strap wedges
    t2, X2, E2 = t, X, None
    Xn = dict(X)
    worst = 1e9
    for k, (n, c) in E.items():
        if n.startswith('PM_avw_bottom') and len(c) == 6:
            cc = [merge[v][0] if v in merge else v for v in c]
            P = [Xn[v] for v in cc]
            vol = min(np.dot(np.cross(P[b] - P[a], P[d] - P[a]), P[e] - P[a]) for a, b, d, e in W3)
            worst = min(worst, vol)
    assert worst > 0, worst
    msg = (f'NOT IN SOURCE (our strap): {len(merge)} strap rim mid-edge nodes merged into the curve\'s refined rim nodes '
           f'(moved {min(d for _, d in merge.values()):.3f}-{max(d for _, d in merge.values()):.3f} mm; '
           + ', '.join(f'{s}->{u}' for s, (u, d) in sorted(merge.items())) + f'); {len(rows) - len(keep)} strap_refine_ties rows removed; '
           f'smallest strap sub-tet {worst:.4f}')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch443.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--tol') + 1]) if '--tol' in a else 0.3)
