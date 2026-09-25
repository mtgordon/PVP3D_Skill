"""Fit each loft's own material to the Abaqus CONN3D2 connectors it replaces (strip model), for every loft.

usage: py -3.10 loft_fit_all.py [MODEL.feb] [--loft NAME ...] [--t 0.49] [--umax 47]
Generalises fit_fan_to_connectors.py (P-arcus only; spacing-based strip widths) to every loft:
  * connector i (length L_i, Abaqus coordinates) becomes a loft strip of length L_i, width w_i and the loft
    thickness t; pulling its far end a distance u gives the uniform stretch 1 + u/L_i (as the connector);
  * w_i: the connectors are ordered along the tissue edge (nearest-neighbour chain through their tissue ends);
    consecutive connectors bound a ruled quad (A_i, B_i, B_j, A_j); connector i owns half of each neighbouring
    quad, w_i = owned area / L_i. All widths are then scaled so that sum w_i L_i = the loft's mesh area (the loft
    is wider than the connector set). For P-arcus this is the spacing method of fit_fan_to_connectors.py;
  * loft pull = sum_i P(1 + u/L_i) w_i t with uniaxial incompressible 1-term Ogden
    P = (c/m)(lam^(m-1) - lam^(-m/2-1)) (FEBio c1 = 2 mu, m1 = alpha); target = sum_i F_i(u), each connector's own
    force-first table (0 in compression, constant past the last point); c by least squares on relative error over
    u in [2, umax] mm (points where the target is < 0.1 % of its maximum are left out), best m on a grid;
  * per behaviour group too (a loft whose connectors have different behaviours: can one material do?);
  * the current loft material (isotropic elastic E 21 = St Venant-Kirchhoff) through the same strips, for scale.
A fan (all connectors from one point, AVW-Para) has extra compliance at the apex that uniform strips ignore.
"""
import argparse
import os
from collections import OrderedDict

import numpy as np

import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR
from loft_survey import survey, tri_area

M_GRID = np.r_[np.arange(1.5, 10.01, 0.25), np.arange(10.5, 40.01, 0.5)]


def ogden_P(lam, c, m):
    return (c / m) * (lam ** (m - 1) - lam ** (-m / 2 - 1))


def svk_P(lam, E):
    return E * lam * (lam ** 2 - 1) / 2


def table_force(tab, u):
    return np.where(u > 0, np.interp(u, tab[:, 1], tab[:, 0]), 0.0)


def order_along_edge(B):
    """Nearest-neighbour chain through the points B, starting at an end of their principal axis."""
    Bc = B - B.mean(axis=0)
    ax = np.linalg.svd(Bc, full_matrices=False)[2][0]
    start = int(np.argmin(Bc @ ax))
    order, left = [start], set(range(len(B))) - {start}
    while left:
        k = min(left, key=lambda j: np.linalg.norm(B[j] - B[order[-1]]))
        order.append(k)
        left.remove(k)
    return order


def quad_area(a, b, c, d):
    return 0.5 * np.linalg.norm(np.cross(c - a, d - b))


def strip_widths(rows, loft_area):
    A = np.array([r['xa'] for r in rows])
    B = np.array([r['xb'] for r in rows])
    L = np.array([r['L'] for r in rows])
    order = order_along_edge(B)
    own = np.zeros(len(rows))
    for i, j in zip(order[:-1], order[1:]):
        q = quad_area(A[i], B[i], B[j], A[j])
        own[i] += q / 2
        own[j] += q / 2
    if len(rows) == 1:
        own[:] = loft_area
    scale = loft_area / own.sum() if own.sum() > 0 else 1.0
    return L, own * scale / L, scale, order


def fit(L, w, t, target, u, wmask):
    lam = 1 + u[:, None] / L[None, :]
    best = None
    for m in M_GRID:
        g = (ogden_P(lam, 1.0, m) * w[None, :] * t).sum(axis=1)
        ok = wmask
        c = np.sum(g[ok] / target[ok]) / np.sum(g[ok] ** 2 / target[ok] ** 2)
        err = np.sqrt(np.mean(((c * g[ok] - target[ok]) / target[ok]) ** 2))
        if best is None or err < best[0]:
            best = (err, c, m, c * g)
    return best


def fit_loft(f, loft, rows, behav, t, umax):
    tabs = [np.array(behav[r['behavior']][0]['table']) for r in rows]
    u = np.linspace(2.0, umax, 91)
    Fi = np.stack([table_force(tb, u) for tb in tabs], axis=1)
    target = Fi.sum(axis=1)
    area = sum(tri_area(f, c) for c in f.elem_blocks[loft][1].values())
    L, w, scale, order = strip_widths(rows, area)
    wmask = target > 1e-3 * target.max()
    err, c, m, pull = fit(L, w, t, target, u, wmask)
    lam = 1 + u[:, None] / L[None, :]
    cur = (svk_P(lam, 21.0) * w[None, :] * t).sum(axis=1)
    print(f"\n== {loft}: {len(rows)} connectors, L {L.min():.1f}-{L.max():.1f} mm, loft area {area:.0f} mm^2 "
          f"(x{scale:.2f} the ruled connector surface), strip widths {w.min():.2f}-{w.max():.2f} mm")
    print(f"   whole loft: c1 {c:.5g} MPa, m1 {m:.2f} (mu0 {c / 2:.4g}, E0 {1.5 * c:.4g} MPa, "
          f"k = 250 mu0 = {125 * c:.4g} MPa), rms rel. error {100 * err:.1f} % over u = 2-{umax:g} mm")
    for x in (2, 5, 10, 20, 30, 40, umax):
        k = int(np.argmin(abs(u - x)))
        print(f"   u = {u[k]:5.1f} mm: connectors {target[k]:10.4g} N | fitted loft {pull[k]:10.4g} N | "
              f"current loft (E 21) {cur[k]:10.4g} N (x{cur[k] / max(target[k], 1e-12):.3g})")
    groups = OrderedDict()
    for i, r in enumerate(rows):
        groups.setdefault(r['behavior'], []).append(i)
    out = dict(loft=loft, n=len(rows), c=c, m=m, err=err, area=area, L=L, w=w, groups={}, order=order)
    if len(groups) > 1:
        for b, idx in groups.items():
            idx = np.array(idx)
            tg = Fi[:, idx].sum(axis=1)
            wm = tg > 1e-3 * tg.max()
            e2, c2, m2, _ = fit(L[idx], w[idx], t, tg, u, wm)
            # stiffness of the group's strips per unit width, relative to the whole-loft fit at u = 10 mm
            k10 = int(np.argmin(abs(u - 10)))
            print(f"   group {b:40s} n {len(idx):2d}: c1 {c2:.4g}, m1 {m2:.2f}, rms {100 * e2:.1f} % | "
                  f"whole-loft material gives {100 * (ogden_P(1 + u[k10] / L[idx], c, m) * w[idx] * t).sum() / max(tg[k10], 1e-12):.0f} % of this group's force at 10 mm")
            out['groups'][b] = (len(idx), c2, m2, e2)
    return out


def fit_all(feb_path, lofts=None, t=0.49, umax=47.0):
    """{loft: fit_loft result} for the named lofts (all lofts when None)."""
    f, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(feb_path)
    by_loft = OrderedDict()
    for fam, (loft, rows) in res.items():
        if loft:
            by_loft.setdefault(loft, []).extend(rows)
    results = OrderedDict()
    for loft, rows in by_loft.items():
        if lofts and loft not in lofts:
            continue
        lnodes = sorted({n for c in f.elem_blocks[loft][1].values() for n in c})
        X = np.array([f.nodes[n] for n in lnodes])
        # a connector belongs to the loft when its tissue end is at a loft node (shared, or coincident when the
        # loft is tied on by tied-node-on-facet, e.g. AVW-Para-L, USL-L), unless a FEBio spring already carries
        # it (PM-middle is the discrete set PM-LA-x%stiff_mat)
        sprung = {frozenset(p) for pairs in f.discsets.values() for p in pairs}
        on = [r for r in rows if np.linalg.norm(X - r['xb'], axis=1).min() < 0.01
              and frozenset((r['fa'], r['fb'])) not in sprung]
        off = [r['elset'] for r in rows if r not in on]
        if off:
            print(f"\n(note: {loft}: connectors not carried by the loft, left out: {off})")
        results[loft] = fit_loft(f, loft, on, behav, t, umax)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('feb', nargs='?', default=os.path.join(RUNS_DIR, 'LPFmsTk_DM1', 'LPFmsTk_DM1.feb'))
    ap.add_argument('--loft', nargs='*')
    ap.add_argument('--t', type=float, default=0.49)
    ap.add_argument('--umax', type=float, default=47.0)
    a = ap.parse_args()
    results = fit_all(a.feb, a.loft, a.t, a.umax)
    print('\n| loft | connectors | behaviours | c1 [MPa] | m1 | rms | k = 250 mu0 |')
    print('|---|---|---|---|---|---|---|')
    for r in results.values():
        print(f"| {r['loft']} | {r['n']} | {len(r['groups']) or 1} | {r['c']:.4g} | {r['m']:.2f} | {100 * r['err']:.1f} % | {125 * r['c']:.4g} |")


if __name__ == '__main__':
    main()
