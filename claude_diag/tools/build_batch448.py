"""Batch 448 (2026-10-23; the user on L446: "there is still a spring that is attached outside of the avw and there are nodes on
the top of the avw that don't have a spring attached"; chose "outer row, every node"). NOT IN SOURCE (as build_batch422 step 3 /
build_batch445): PM_conn (the AVW -> top of the PM arch) gets one spring on EVERY node of the AVW's outer distal row (the row
its springs already use, seam nodes with the curves included: 27 nodes), paired in x order with the fixed PM ends picked evenly
by index (26 ends for 27 springs: one end takes 2). PM_conn_mat <scale> x N_old / N_new so the family's total force is unchanged.
The middle / inner rows of the end face stay free. Old patch anchors stay (no spring on them).
usage: py -3.10 build_batch448.py SRC_RUN NEW_RUN"""
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
    pairs = [(int(p), int(q)) for p, q in re.findall(r'<delem>(\d+),(\d+)</delem>', m.group(2))]
    assert len(pairs) == len(re.findall('<delem', m.group(2)))
    near = [p if q in fixed else q for p, q in pairs]
    far = sorted((q if q in fixed else p for p, q in pairs), key=lambda v: X[v][0])
    wend = [v for v in near if dom.get(v) == {'_PickedSet347'}]
    y0, z0 = np.median([X[v][1] for v in wend]), np.median([X[v][2] for v in wend])
    row = sorted((v for v in X if '_PickedSet347' in dom.get(v, ()) and abs(X[v][1] - y0) < 0.06 and abs(X[v][2] - z0) < 0.12),
                 key=lambda v: X[v][0])
    fpick = [far[i] for i in np.round(np.linspace(0, len(far) - 1, len(row))).astype(int)]
    ind = re.match(r'\s*', m.group(2)).group(0)
    body = ''.join(f'{ind}<delem>{n},{f}</delem>' for n, f in zip(row, fpick)) + re.search(r'\s*$', m.group(2)).group(0)
    t = t[:m.start(2)] + body + t[m.end(2):]
    dm = re.search(r'<discrete dmat="(\d+)" discrete_set="PM_conn" />', t).group(1)
    assert len(re.findall(r'<discrete dmat="%s" ' % dm, t)) == 1, 'PM_conn_mat shared'
    mats = re.findall(r'<discrete_material id="\d+"[^>]*>.*?</discrete_material>', t, re.S)
    mat = mats[int(dm) - 1]
    assert 'name="PM_conn_mat"' in mat, mat[:80]
    s0 = float(re.search(r'<scale>([^<]+)</scale>', mat).group(1)); s1 = s0 * len(pairs) / len(row)
    t = t.replace(mat, re.sub(r'<scale>[^<]+</scale>', f'<scale>{s1:.6g}</scale>', mat))
    off = [v for v in near if dom.get(v) != {'_PickedSet347'}]
    dup = sorted({v for v in fpick if fpick.count(v) > 1})
    msg = (f'NOT IN SOURCE (spring ends): PM_conn {len(pairs)} -> {len(row)} springs, one on every node of the AVW outer distal row '
           f'(x {X[row[0]][0]:.1f}..{X[row[-1]][0]:.1f}, seam nodes included), paired in x order with the fixed PM ends (end {dup} takes 2); '
           f'PM_conn_mat scale {s0:g} -> {s1:.6g} (family total unchanged); {len(off)} ends moved off the AVW / patches {off}; '
           f'middle / inner rows of the end face free')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch448.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
