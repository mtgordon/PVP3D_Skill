"""Batch 405 (2026-10-22; the user: a fillet at the curve's distal end, "do both so I can see them"). Node moves only (no
connectivity change; ties tie displacements, so they keep working), NOT IN SOURCE (the curves are ours):
  --mode rim   : round the END-FACE RIM (the edge between the curves' distal end face and their sides; it folds from 71-94 deg
                 to 143-168 deg by t 1 in L393): every rim node (a node of both an end face and a side face of curve_all, end
                 faces = within 1.5 mm of the end along the end tangent with their normal along it) moves into the material
                 along the bisector -(n_end + n_side) by DELTA (= 0.414 R, the corner of a radius-R fillet). Rim nodes shared
                 with the strap / AVW move too (the user's yes).
  --mode ridge : round the LENGTHWISE creases (curve_all edges along the curve, dihedral > ANG) over the last ZLEN mm: every
                 free (curve-only) node on such a crease moves ALPHA x w(s) of the way to the midpoint of its two neighbours
                 around the cross-section (w = 1 to ZLEN - 2 mm from the end, ramping to 0 at ZLEN).
Each move is cut back (bisection) until every element containing the node keeps >= QMIN of its rest volume.
usage: py -3.10 build_batch405.py SRC_RUN NEW_RUN --mode rim|ridge [--r 1.0] [--zlen 5] [--alpha 0.5] [--ang 25] [--qmin 0.6] [--dry]"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402

CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]
W3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def vol(X, c):
    sp = T5 if len(c) == 8 else (W3 if len(c) == 6 else None)
    if sp is None:
        return None
    return sum(np.dot(np.cross(X[c[b]] - X[c[a]], X[c[d]] - X[c[a]]), X[c[e]] - X[c[a]]) / 6 for a, b, d, e in sp)


def setup(t, X, E):
    a = t.index('<Surface name="curve_all">'); b = t.index('</Surface>', a)
    F = [[int(v) for v in s.split(',')] for s in re.findall(r'>([0-9,]+)</', t[a:b])]
    F = [f for f in F if len(f) == 4]
    nrm = []
    for f in F:
        P = np.array([X[v] for v in f]); n = np.cross(P[2] - P[0], P[3] - P[1]); nrm.append(n / np.linalg.norm(n))
    side = {}
    for sg in (-1, 1):
        cn = np.array([X[v] for nm, c in E.values() if nm in CURVES for v in c if np.sign(X[v][0]) == sg])
        lo = cn[:, 2].min(); e = cn[cn[:, 2] < lo + 3].mean(0); u = cn[(cn[:, 2] > lo + 8) & (cn[:, 2] < lo + 12)].mean(0)
        tan = (e - u) / np.linalg.norm(e - u)
        smax = max(X[v] @ tan for nm, c in E.values() if nm in CURVES for v in c if np.sign(X[v][0]) == sg)
        side[sg] = (tan, smax)
    return F, nrm, side


def main(src, new, mode, R, zlen, alpha, ang, qmin, dry):
    t, X0, E = load(src)
    X = {k: v.copy() for k, v in X0.items()}
    F, nrm, side = setup(t, X0, E)
    owner = defaultdict(list)
    for k, (nm, c) in E.items():
        for v in c:
            owner[v].append(k)
    dn = defaultdict(set)
    for nm, c in E.values():
        for v in c:
            dn[v].add(nm)
    target = {}
    if mode == 'rim':
        endf, sidef = defaultdict(list), defaultdict(list)
        for f, n in zip(F, nrm):
            sg = np.sign(np.mean([X0[v][0] for v in f])); tan, smax = side[sg]
            if all(smax - X0[v] @ tan < 1.5 for v in f) and n @ tan > 0.8:
                for v in f:
                    endf[v].append(n)
            else:
                for v in f:
                    sidef[v].append(n)
        for v in set(endf) & set(sidef):
            ne = np.mean(endf[v], 0); ns = np.mean(sidef[v], 0)
            ne /= np.linalg.norm(ne); ns /= np.linalg.norm(ns)
            d = -(ne + ns); d /= np.linalg.norm(d)
            target[v] = X0[v] + 0.414 * R * d
    else:
        edges = defaultdict(list)
        for i, f in enumerate(F):
            for k in range(4):
                edges[frozenset((f[k], f[(k + 1) % 4]))].append(i)
        nb = defaultdict(set)   # around-the-section surface neighbours
        crease = set()
        for e_, fs in edges.items():
            p, q = tuple(e_)
            sg = np.sign(X0[p][0]); tan, smax = side[sg]
            dvec = X0[q] - X0[p]; along = abs(dvec @ tan) / np.linalg.norm(dvec)
            if along < 0.4:
                nb[p].add(q); nb[q].add(p)
            if len(fs) == 2 and along > 0.7:
                da = np.degrees(np.arccos(np.clip(nrm[fs[0]] @ nrm[fs[1]], -1, 1)))
                if da > ang and smax - max(X0[p] @ tan, X0[q] @ tan) < zlen:
                    crease |= {p, q}
        for v in crease:
            if dn[v] - set(CURVES) or len(nb[v]) != 2:
                continue
            sg = np.sign(X0[v][0]); tan, smax = side[sg]
            s = smax - X0[v] @ tan
            w = 1.0 if s < zlen - 2 else max(0.0, (zlen - s) / 2)
            m = np.mean([X0[q] for q in nb[v]], 0)
            target[v] = X0[v] + alpha * w * (m - X0[v])
    # cut back each move until every element containing the node keeps >= qmin of its rest volume
    rest = {k: vol(X0, c) for k, (nm, c) in E.items()}
    moved, cut = {}, 0
    for v in sorted(target, key=lambda v: -np.linalg.norm(target[v] - X0[v])):
        lo_, hi_ = 0.0, 1.0
        def ok(s_):
            X[v] = X0[v] + s_ * (target[v] - X0[v])
            return all(rest[k] is None or vol(X, E[k][1]) >= qmin * rest[k] for k in owner[v])
        if ok(1.0):
            s_ = 1.0
        else:
            for _ in range(20):
                mid = (lo_ + hi_) / 2
                lo_, hi_ = (mid, hi_) if ok(mid) else (lo_, mid)
            s_ = lo_; cut += 1
        X[v] = X0[v] + s_ * (target[v] - X0[v])
        if s_ > 0:
            moved[v] = np.linalg.norm(X[v] - X0[v])
    mv = np.array(list(moved.values())) if moved else np.zeros(1)
    shared = sum(1 for v in moved if dn[v] - set(CURVES))
    print(f'{mode}: {len(target)} nodes targeted, {len(moved)} moved (median {np.median(mv):.2f}, max {mv.max():.2f} mm), '
          f'{cut} cut back for element quality (>= {qmin:g} of rest volume), {shared} of them shared with another domain')
    if dry:
        return X0, X, moved
    def node_sub(m):
        k = int(m.group(1))
        return f'{m.group(0)[:m.start(2) - m.start(0)]}{X[k][0]:.7g},{X[k][1]:.7g},{X[k][2]:.7g}</node>' if k in moved else m.group(0)
    t = re.sub(r'<node id="(\d+)">([^<]+)</node>', node_sub, t)
    what = (f'the curves\' distal END-FACE RIM rounded: {len(moved)} rim nodes moved into the material along the end / side '
            f'bisector (target 0.414 x R {R:g} mm)' if mode == 'rim' else
            f'the LENGTHWISE creases (dihedral > {ang:g} deg) over the last {zlen:g} mm rounded: {len(moved)} free crease nodes moved '
            f'{alpha:g} x w(s) toward their cross-section neighbours\' midpoint')
    msg = (f'NOT IN SOURCE (the user 10-22, fillet): {what}; moves median {np.median(mv):.2f} / max {mv.max():.2f} mm, {cut} cut back '
           f'to keep every element >= {qmin:g} of its rest volume; {shared} moved nodes shared with the strap / AVW; node '
           f'positions only (connectivity, ties, surfaces unchanged)')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch405.py)\n{msg}\n')
    return X0, X, moved


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], a[a.index('--mode') + 1], g('--r', 1.0), g('--zlen', 5), g('--alpha', 0.5), g('--ang', 25), g('--qmin', 0.6),
         '--dry' in a)
