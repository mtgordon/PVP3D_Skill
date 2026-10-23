"""Batch 445 (2026-10-23; the user: the springs from the top of the PM arch "just attach to the avw for now (not the corners)").
NOT IN SOURCE (as build_batch422 step 3): PM_conn (26 springs, the AVW -> top of the PM arch) had 4 ends off the AVW interior:
2 on spring patches at the AVW-curve seam nodes of the distal row and 2 on the curves' outer top corners (patches mostly on the
strap sheet, L440g on). Here all 26 AVW-side ends go onto the AVW's own distal-row nodes, seam nodes excluded (25 nodes, picked
evenly by index, so one node takes 2 springs), paired in x order with the fixed PM ends (no crossing). Spring material, far ends
unchanged; the old patch anchors stay (no spring on them).
usage: py -3.10 build_batch445.py SRC_RUN NEW_RUN"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new):
    t, X, E = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    dom = {}
    for n, c in E.values():
        for v in c:
            dom.setdefault(v, set()).add(n)
    fixed = set()
    for mm in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', t):
        if re.search(r'<bc [^>]*node_set="%s"[^>]*type="zero displacement"' % re.escape(mm.group(1)), t):
            fixed |= {int(v) for v in re.findall(r'\d+', mm.group(2))}
    m = re.search(r'(<DiscreteSet name="PM_conn">)(.*?)(</DiscreteSet>)', t, re.S)
    rows = re.findall(r'(\s*<delem[^>]*>)(\d+),(\d+)(</delem>)', m.group(2))
    pairs = [(int(p), int(q)) for _, p, q, _ in rows]
    near = [p if q in fixed else q for p, q in pairs]
    far = sorted((q if q in fixed else p for p, q in pairs), key=lambda v: X[v][0])
    wend = [v for v in near if dom.get(v) == {'_PickedSet347'}]
    y0, z0 = np.median([X[v][1] for v in wend]), np.median([X[v][2] for v in wend])
    row = sorted((v for v in X if '_PickedSet347' in dom.get(v, ()) and abs(X[v][1] - y0) < 0.06 and abs(X[v][2] - z0) < 0.12),
                 key=lambda v: X[v][0])
    inner = [v for v in row if dom[v] == {'_PickedSet347'}]
    seam = [v for v in row if v not in inner]
    assert len(seam) == 2 and seam[0] == row[0] and seam[1] == row[-1], (seam, row[0], row[-1])
    pick = [inner[i] for i in np.round(np.linspace(0, len(inner) - 1, len(far))).astype(int)]
    off = [v for v in near if v not in inner]
    body = ''.join(f'{o}{n},{f}{c}' for (o, _, _, c), n, f in zip(rows, pick, far)) + re.search(r'\s*$', m.group(2)).group(0)
    t = t[:m.start(2)] + body + t[m.end(2):]
    dup = sorted({v for v in pick if pick.count(v) > 1})
    msg = (f'NOT IN SOURCE (spring ends): PM_conn {len(far)} AVW-side ends onto the AVW distal row only ({len(inner)} interior nodes '
           f'x {X[inner[0]][0]:.1f}..{X[inner[-1]][0]:.1f}, seam nodes {seam} excluded, node {dup} takes 2), paired in x order with the '
           f'fixed PM ends; {len(off)} ends moved off the seam / corner patches {off}')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch445.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
