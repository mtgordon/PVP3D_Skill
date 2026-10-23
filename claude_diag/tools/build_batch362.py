"""Batch 362 (2026-10-20 night; the user's go (2) of NEXT_SESSION_PROMPT_2026-10-20b: "try the new geometry" = option (2) of
curve_end_view/curve_end_proposals.png, the distal curve pre-bent 4 mm outward). In L352 the distal ~25 mm of each curve (rest z -40
to -15) is squeezed lengthwise between the descending wall and the dropping perineal body and buckles outward (+y) into a hairpin.
NOT IN SOURCE (rest geometry changed): every node of the canal domains (DOMS: AVW / PVW hexes and the curves) at rest z in
(Z0, Z1) is moved by u_y = AMP * sin^2(pi (z - Z0) / (Z1 - Z0)) (+y on both sides: the two curves are mirror images about x = 0, so
+y is "outward" for both, the way they buckle). The curve's cross-section spans from the AVW's lateral edge to the PVW's (it is the
lateral sulcus), so the curve cannot move without the walls it joins: the whole distal canal slab moves by the same profile, a
smooth bend of the canal end with both ends unchanged. sin^2 (the clamped-clamped first buckling mode) rather than sin: no corner
at the wall join above or the perineal body below (sin would tilt the ends ~27 deg). Nodes of other parts move only when shared
with a canal element (the strap's curve-edge nodes, near z -38 to -40, where u_y < 0.2 mm). Every linear constraint touching a
moved node keeps its weights; its dependent node (weight 1) is put back at the weighted master position plus its rest residual,
so every tie is exactly as it was (iterated for chains).
Checks printed: moved nodes per domain, max shift, tie residuals, hex8 / penta6 corner Jacobians of every element touching a moved
node (min scaled Jacobian before / after), fixed-BC nodes moved (must be 0), and the walls <-> LA rest distance (offset 2 mm).
usage: py -3.10 build_batch362.py SRC_RUN NEW_RUN [--amp 4] [--z0 -40] [--z1 -15] [--dry]"""
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from build_batch162 import write  # noqa: E402
from build_batch361 import elements  # noqa: E402
from feb_model import pt_tri_dist  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
DOMS = ('_PickedSet347', '_PickedSet64', 'canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
# corner -> its three edge neighbours (right-handed for a positive element)
HEXC = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7), (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]
PENC = [(0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 1, 5), (3, 5, 4, 0), (4, 3, 5, 1), (5, 4, 3, 2)]


def min_sj(P, ty):
    out = 9.0
    for c, a, b, d in (HEXC if ty == 'hex8' else PENC):
        u, v, w = P[a] - P[c], P[b] - P[c], P[d] - P[c]
        out = min(out, np.dot(np.cross(u, v), w) / (np.linalg.norm(u) * np.linalg.norm(v) * np.linalg.norm(w) + 1e-30))
    return out


def constraints(t):
    out = []   # (block name, start, end, [(node, weight)] of the x row)
    for blk in re.finditer(r'<constraint name="([^"]+)" type="linear constraint">(.*?)</constraint>', t, re.S):
        for lc in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', blk.group(2), re.S):
            rows = re.findall(r'<node id="(\d+)" bc="(\w)">([^<]+)</node>', lc)
            if rows[0][1] == 'x':
                out.append((blk.group(1), [(int(n), float(w)) for n, b, w in rows]))
    return out


def main(src, new, amp, z0, z1, dry):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    E = elements(t)
    canal = set(v for d in DOMS for e in E[d][1] for v in e)
    Y = {k: v.copy() for k, v in X.items()}
    moved = set()
    for v in canal:
        z = X[v][2]
        if z0 < z < z1:
            Y[v][1] += amp * np.sin(np.pi * (z - z0) / (z1 - z0)) ** 2
            moved.add(v)
    # ties: dependent = weighted masters + rest residual
    C = constraints(t)
    res0 = []
    for nm, rows in C:
        (d, w1), ms = rows[0], rows[1:]
        assert abs(w1 - 1) < 1e-9
        res0.append(X[d] + sum(w * X[n] for n, w in ms))   # weights of the masters are negative
    tied_dep = set()
    for it in range(10):
        ch = 0
        for (nm, rows), r0 in zip(C, res0):
            (d, _), ms = rows[0], rows[1:]
            if not (any(n in moved for n, _ in ms) or d in moved):
                continue
            p = r0 - sum(w * Y[n] for n, w in ms)
            if np.linalg.norm(p - Y[d]) > 1e-9:
                ch += 1; Y[d] = p; moved.add(d); tied_dep.add(d)
        if not ch:
            break
    print(f'ties settled after {it + 1} passes')
    res = max(np.linalg.norm((Y[rows[0][0]] + sum(w * Y[n] for n, w in rows[1:])) - r0) for (nm, rows), r0 in zip(C, res0))
    shift = {v: np.linalg.norm(Y[v] - X[v]) for v in moved}
    dev = max((abs(np.linalg.norm(Y[v] - X[v]) - amp * np.sin(np.pi * (X[v][2] - z0) / (z1 - z0)) ** 2) if z0 < X[v][2] < z1 else
               np.linalg.norm(Y[v] - X[v])) for v in tied_dep) if tied_dep else 0
    print(f'moved nodes {len(moved)} (canal {len(moved & canal)}, other {len(moved - canal)}); max shift {max(shift.values()):.3f} mm; '
          f'tie dependents re-placed {len(tied_dep)} (max off the profile {dev:.3f} mm); max tie residual change {res:.2e} mm')
    from collections import Counter
    print('moved by domain:', dict(Counter(d for d, (ty, el) in E.items() for v in set(x for e in el for x in e) if v in moved)))
    # element quality
    worst = []
    for d, (ty, el) in E.items():
        if ty not in ('hex8', 'penta6'):
            continue
        for e in el:
            if not any(v in moved for v in e):
                continue
            a = min_sj(np.array([X[v] for v in e]), ty); b = min_sj(np.array([Y[v] for v in e]), ty)
            worst.append((b, a, d, e))
    worst.sort(key=lambda r: r[0])
    print(f'solid elements touching moved nodes {len(worst)}; min scaled corner Jacobian after {worst[0][0]:.3f} (before {worst[0][1]:.3f}, '
          f'{worst[0][2]}); negative after {sum(1 for r in worst if r[0] <= 0)} (before {sum(1 for r in worst if r[1] <= 0)})')
    by = {}
    for b, a, d, e in worst:
        by.setdefault(d, []).append((b, a))
    for d, v in by.items():
        print(f'   {d}: {len(v)} elements, min sJ after {min(x[0] for x in v):.3f} / before {min(x[1] for x in v):.3f}; '
              f'largest drop {max(x[1] - x[0] for x in v):.3f}')
    # fixed BCs
    fixed = set()
    for ns in re.findall(r'<bc name="[^"]+" node_set="([^"]+)" type="zero', t):
        m = re.search(rf'<NodeSet name="{re.escape(ns)}">([^<]*)</NodeSet>', t)
        fixed |= {int(v) for v in re.findall(r'\d+', m.group(1))}
    print(f'fixed-BC nodes moved: {len(fixed & moved)}')
    # walls <-> LA rest distance (sliding-elastic offset 2 mm, search radius 3)
    for prim, sec in (('walls_LA_primary', 'PVW_LA_secondary'), ('PVW_LA_primary', 'PVW_LA_secondary')):
        fp = re.search(rf'<Surface name="{prim}">(.*?)</Surface>', t, re.S).group(1)
        fs = re.search(rf'<Surface name="{sec}">(.*?)</Surface>', t, re.S).group(1)
        pn = sorted(set(int(v) for c in re.findall(r'>([0-9,]+)</', fp) for v in c.split(',')) & moved)
        tris = []
        for c in re.findall(r'>([0-9,]+)</', fs):
            q = [int(v) for v in c.split(',')]
            tris += [(q[0], q[1], q[2])] + ([(q[0], q[2], q[3])] if len(q) == 4 else [])
        for lab, P in (('before', X), ('after', Y)):
            tc = cKDTree(np.array([np.mean([P[v] for v in tr], 0) for tr in tris]))
            dm = []
            for v in pn:
                dm.append(min(pt_tri_dist(P[v], *(P[u] for u in tris[j])) for j in tc.query(P[v], 12)[1]))
            dm = np.array(dm)
            print(f'   {prim} moved nodes {len(pn)} -> {sec}: {lab} min {dm.min():.2f} mm, < 2 mm (offset) {int((dm < 2).sum())}')
    if dry:
        return
    def sub(m):
        i = int(m.group(1))
        if i in moved:
            return f'<node id="{i}">{Y[i][0]:.6f},{Y[i][1]:.6f},{Y[i][2]:.6f}</node>'
        return m.group(0)
    t2 = re.sub(r'<node id="(\d+)">[^<]+</node>', sub, t)
    msg = (f'NOT IN SOURCE (rest geometry): the distal canal pre-bent {amp:g} mm in +y (option 2 of curve_end_proposals.png): every '
           f'node of {", ".join(DOMS)} at rest z in ({z0:g}, {z1:g}) moved by u_y = {amp:g} sin^2(pi (z - {z0:g}) / {z1 - z0:g}) '
           f'({len(moved & canal)} nodes); {len(tied_dep)} tie-dependent nodes re-placed on their masters (weights unchanged, '
           f'residuals as at rest); {len(moved - canal)} non-canal nodes moved (shared / tied); min scaled corner Jacobian of the '
           f'touched solids {worst[0][0]:.3f} (was {worst[0][1]:.3f}); fixed-BC nodes moved {len(fixed & moved)} (tools/build_batch362.py)')
    print(msg)
    write(new, t2, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch362.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: float(a[a.index(k_) + 1]) if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--amp', 4.0), g('--z0', -40.0), g('--z1', -15.0), '--dry' in a)
