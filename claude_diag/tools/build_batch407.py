"""Batch 407 (2026-10-22; the user: "I also would like to see a thicker strap"). The fan strap (PM_avw_bottom_left/right_sheet_solid,
one penta6 layer ~0.49 mm, extruded from its base sheet back INTO the curve) is flipped to the distal side of its base (as
build_batch403's step 1, no cap) at thickness TH: top' = base - TH x unit(top - base); every wedge renumbered (c3, c4, c5, c0, c1, c2)
so it stays positive; all-strap surface facets re-wound. To keep the membrane stiffness (thickness x modulus) the two strap
materials' matrix c1 and fibre E are scaled by TH0 / TH (TH0 = the median rest thickness). A top node that would sit inside
another solid at rest (next to the PeB) is pulled back toward its base (0.05 mm margin, at least 0.2 of TH). The base nodes
(rim tie, junction / patch ties, PM-end anchors), element ids, fibre directions are unchanged. NOT IN SOURCE (the strap is ours).
Note: the bending stiffness still rises ~(TH / TH0)^2.
usage: py -3.10 build_batch407.py SRC_RUN NEW_RUN [--th 1.0]"""
import os
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from build_batch405 import vol, T5, W3  # noqa: E402
from junction_view_data import load  # noqa: E402

STRAP = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
MATS = ('PM_avw_bottom_left_sheet_mat', 'PM_avw_bottom_right_sheet_mat')


def main(src, new, th):
    t, X0, E = load(src)
    X = {k: v.copy() for k, v in X0.items()}
    dn = defaultdict(set)
    for nm, c in E.values():
        for v in c:
            dn[v].add(nm)
    wed = {k: c for k, (nm, c) in E.items() if nm in STRAP}
    partner = {}
    for c in wed.values():
        for i in range(3):
            assert partner.get(c[i], c[i + 3]) == c[i + 3]
            partner[c[i]] = c[i + 3]
    top = set(partner.values())
    assert all(dn[v] <= set(STRAP) for v in top)
    th0 = float(np.median([np.linalg.norm(X0[t_] - X0[b_]) for b_, t_ in partner.items()]))
    for b_, t_ in partner.items():
        u = X0[t_] - X0[b_]; X[t_] = X0[b_] - th * u / np.linalg.norm(u)
    newwed = {k: [c[3], c[4], c[5], c[0], c[1], c[2]] for k, c in wed.items()}
    # keep the top nodes outside other solids
    tets = []
    for k, (nm, c) in E.items():
        if nm in STRAP:
            continue
        sp = T5 if len(c) == 8 else (W3 if len(c) == 6 else None)
        if sp:
            tets += [[c[i] for i in s] for s in sp]
    TP = np.array([[X0[v] for v in q] for q in tets]); cen = TP.mean(1)
    tr = cKDTree(cen); RM = np.max(np.linalg.norm(TP - cen[:, None], axis=2))

    def inside(p, skip):
        for j in tr.query_ball_point(p, RM):
            if skip & set(tets[j]):
                continue
            a_, b_, c_, d_ = TP[j]
            lam = np.linalg.solve(np.column_stack([b_ - a_, c_ - a_, d_ - a_]), p - a_)
            if min(1 - lam.sum(), lam.min()) > 1e-6:
                return True
        return False
    short = []
    for b_, t_ in partner.items():
        # (L408 check) sample the whole extrusion, not only its end: a thick strap can pass THROUGH a thin lip of the PeB
        off = X[t_] - X0[b_]
        hits = [f for f in np.linspace(0.05, 1, 20) if inside(X0[b_] + f * off, {b_, t_})]
        if not hits:
            continue
        s_ = max(hits[0] - 0.05 - 0.05 / np.linalg.norm(off), 0.2)
        X[t_] = X0[b_] + s_ * off; short.append(round(s_, 2))
    bad = [k for k, c in newwed.items() if vol(X, c) <= 0]
    assert not bad, f'non-positive wedges {bad}'
    def node_sub(m):
        k = int(m.group(1))
        return f'{m.group(0)[:m.start(2) - m.start(0)]}{X[k][0]:.7g},{X[k][1]:.7g},{X[k][2]:.7g}</node>' if k in top else m.group(0)
    t = re.sub(r'<node id="(\d+)">([^<]+)</node>', node_sub, t)

    def el_sub(m):
        k = int(m.group(1))
        return f'{m.group(0)[:m.start(2) - m.start(0)]}{",".join(map(str, newwed[k]))}</elem>\n' if k in newwed else m.group(0)
    t = re.sub(r'\t\t\t<elem id="(\d+)">([0-9,]+)</elem>\n', el_sub, t)
    scale = th0 / th
    for mname in MATS:
        a = t.index(f'name="{mname}"'); b = t.index('</material>', a)
        blk = t[a:b]
        blk = re.sub(r'<c1>([^<]+)</c1>', lambda m: f'<c1>{float(m.group(1)) * scale:.6g}</c1>', blk)
        blk = re.sub(r'<E>([^<]+)</E>', lambda m: f'<E>{float(m.group(1)) * scale:.6g}</E>', blk)
        t = t[:a] + blk + t[b:]
    strapn = set(partner) | top
    nrew = defaultdict(int)

    def surf_sub(ms):
        name, body = ms.group(1), ms.group(2)
        out = []
        for g, fid, nodes in re.findall(r'\t\t\t<(quad4|tri3) id="(\d+)">([0-9,]+)</\1>\n', body):
            q = [int(v) for v in nodes.split(',')]
            if set(q) <= strapn:
                q = q[::-1]; nrew[name] += 1
            out.append((g, q))
        return (f'\t\t<Surface name="{name}">\n' + ''.join(f'\t\t\t<{g} id="{i}">{",".join(map(str, q))}</{g}>\n'
                                                       for i, (g, q) in enumerate(out, 1)) + '\t\t</Surface>\n')
    t = re.sub(r'\t\t<Surface name="([^"]+)">\n(.*?)\t\t</Surface>\n', surf_sub, t, flags=re.S)
    msg = (f'NOT IN SOURCE (the user 10-22): the fan strap flipped to the distal side of its base at thickness {th:g} mm (rest '
           f'{th0:.2f}); matrix c1 and fibre E x{scale:.3f} (membrane stiffness kept, bending ~x{(th / th0) ** 2:.1f}); '
           f'{len(wed)} wedges renumbered; {len(short)} top nodes next to other solids pulled back (to {min(short or [1]):.2f}-1 of TH); '
           f'facets re-wound: ' + ', '.join(f'{k} {v}' for k, v in nrew.items()))
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch407.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--th') + 1]) if '--th' in a else 1.0)
