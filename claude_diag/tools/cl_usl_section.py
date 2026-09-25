"""CL/USL loft material from the Abaqus connectors by cross-section (user, 2026-09-24): "make the CL/USL tissue properties
similar to Abaqus, accounting for the number of connectors it had, compared with the cross-sectional area it has now".

usage: py -3.10 cl_usl_section.py [MODEL.feb] [--t 0.49]
Per loft (CL-L/R, USL-L/R):
  * the connectors it replaces: count per behaviour group, lengths; total force F_tot(u) = sum of every connector's own
    force-first table at elongation u (0 at u <= 0, constant past the last point);
  * the loft's cross-section: width x thickness t, the width measured three ways: mesh area / mean connector length
    (the average width, used for the property), the anchor-edge length (the CL/USL line side) and the tissue-edge
    length (the vaginal-wall side), each along the loft's boundary loop;
  * the tissue curve: nominal stress sigma = F_tot(u) / A (A from the average width) at nominal strain eps = u / L_mean;
  * a 1-term and a 2-term Ogden fit to sigma(eps) (FEBio c_i, m_i; uniaxial, incompressible), compared with the strip-
    model fit (loft_fit_all.py) at the same strains.
Writes claude_diag/cl_usl_section_<loft>.txt as Abaqus-style *Uniaxial Test Data (nominal stress, nominal strain) for
scripts/fit_yeoh.py / fit_ogden.pl as well.
"""
import argparse
import os
from collections import Counter, defaultdict

import numpy as np

import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR, DIAG_DIR
from loft_survey import survey, tri_area
from loft_fit_all import fit_all, ogden_P

FAMS = ('CL-L', 'CL-R', 'USL-L', 'USL-R')


def boundary_loop(conns):
    cnt = defaultdict(int)
    for c in conns:
        for k in range(len(c)):
            cnt[tuple(sorted((c[k], c[(k + 1) % len(c)])))] += 1
    adj = defaultdict(list)
    for (a, b), n in cnt.items():
        if n == 1:
            adj[a].append(b)
            adj[b].append(a)
    s = min(adj)
    loop, prev, cur = [s], None, s
    while True:
        nxt = [n for n in adj[cur] if n != prev and n != s]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        loop.append(cur)
    return loop


def arc_length(loop, members, X):
    """Length of the boundary arc that spans all `members` (the complement of the largest gap between them)."""
    I = [i for i, n in enumerate(loop) if n in members]
    if len(I) < 2:
        return 0.0, 0
    gaps = [((I[(k + 1) % len(I)] - I[k]) % len(loop), k) for k in range(len(I))]
    big = max(gaps)[1]
    start = I[(big + 1) % len(I)]
    arc = [loop[(start + j) % len(loop)] for j in range((I[big] - start) % len(loop) + 1)]
    return float(sum(np.linalg.norm(X[b] - X[a]) for a, b in zip(arc[:-1], arc[1:]))), len(I)


def fit_ogden(eps, sig, nterm):
    """Least squares on relative error; 1 term: c by closed form on an m grid; 2 terms: (m1, m2) grid, c by NNLS-like
    2x2 solve with c >= 0."""
    lam = 1 + eps
    ok = sig > 1e-3 * sig.max()
    grid = np.r_[np.arange(1.5, 10.01, 0.25), np.arange(10.5, 40.01, 0.5)]
    best = None
    if nterm == 1:
        for m in grid:
            g = ogden_P(lam, 1.0, m)
            c = np.sum(g[ok] / sig[ok]) / np.sum(g[ok] ** 2 / sig[ok] ** 2)
            err = np.sqrt(np.mean(((c * g[ok] - sig[ok]) / sig[ok]) ** 2))
            if best is None or err < best[0]:
                best = (err, [(c, m)])
        return best
    for i, m1 in enumerate(grid):
        for m2 in grid[i + 1:]:
            G = np.stack([ogden_P(lam, 1.0, m1), ogden_P(lam, 1.0, m2)], axis=1)[ok] / sig[ok][:, None]
            c, *_ = np.linalg.lstsq(G, np.ones(ok.sum()), rcond=None)
            if (c < 0).any():
                continue
            err = np.sqrt(np.mean((G @ c - 1) ** 2))
            if best is None or err < best[0]:
                best = (err, [(c[0], m1), (c[1], m2)])
    return best


def section_fits(feb_path, t=0.49, strip=False):
    """{loft: dict} for the CL/USL lofts: the section curve (u, F_tot, eps, sig), the geometry (Lm, area, widths, A) and
    the 1- and 2-term Ogden fits (e1, t1, e2, t2: rms error and [(c, m)]); with strip=True also the strip-model fit."""
    f, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(feb_path)
    strips = fit_all(feb_path, [res[fam][0] for fam in FAMS]) if strip else {}
    out = {}
    for fam in FAMS:
        loft, rows = res[fam]
        tabs = [np.array(behav[r['behavior']][0]['table']) for r in rows]
        L = np.array([r['L'] for r in rows])
        Lm = float(L.mean())
        conn = list(f.elem_blocks[loft][1].values())
        area = float(sum(tri_area(f, c) for c in conn))
        loop = boundary_loop(conn)
        anchor, _ = arc_length(loop, {n for n in loop if n in bcs}, f.nodes)
        tissue, _ = arc_length(loop, {n for n in loop if elem_of[n] - {loft}}, f.nodes)
        w = area / Lm
        A = w * t
        umax = max(tb[-1, 1] for tb in tabs)
        u = np.linspace(0.5, umax, 200)
        Ft = sum(np.where(u > 0, np.interp(u, tb[:, 1], tb[:, 0]), 0.0) for tb in tabs)
        eps, sig = u / Lm, Ft / A
        e1, t1 = fit_ogden(eps, sig, 1)
        e2, t2 = fit_ogden(eps, sig, 2)
        out[loft] = dict(fam=fam, rows=rows, grp=Counter(r['behavior'] for r in rows), behav=behav, L=L, Lm=Lm,
                         area=area, anchor=anchor, tissue=tissue, w=w, A=A, u=u, Ft=Ft, eps=eps, sig=sig, e1=e1, t1=t1,
                         e2=e2, t2=t2, strip=strips.get(loft))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('feb', nargs='?', default=os.path.join(RUNS_DIR, 'L6_ls', 'L6_ls.feb'))
    ap.add_argument('--t', type=float, default=0.49)
    a = ap.parse_args()
    out = section_fits(a.feb, a.t, strip=True)
    print('\n| loft | connectors (groups) | L min / mean / max [mm] | loft area [mm^2] | width: average / anchor edge / tissue '
          'edge [mm] | cross-section A (average) [mm^2] | F_tot at u = 25 / 36.5 / 51.9 mm [N] | stress at those u [MPa] '
          '(strain) | Ogden 1-term c1, m1 (rms) | Ogden 2-term | strip fit c1, m1 |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for loft, r in out.items():
        rows, grp, behav, L, Lm, area, anchor, tissue, w, A, u, Ft, eps, sig, e1, t1, e2, t2, s = (
            r[k] for k in ('rows', 'grp', 'behav', 'L', 'Lm', 'area', 'anchor', 'tissue', 'w', 'A', 'u', 'Ft', 'eps',
                           'sig', 'e1', 't1', 'e2', 't2', 'strip'))
        umax = float(u[-1])
        pts = [25.0, 36.5, umax]
        Fp = [float(np.interp(p, u, Ft)) for p in pts]
        sp = [float(np.interp(p, u, sig)) for p in pts]
        f1 = {b: float(behav[b][0]['table'][1][0]) for b in grp}          # each group's force at its first point
        top = max(f1.values())
        glabel = ', '.join(f'{v} x{f1[b] / top:.2g}' for b, v in grp.items())
        print(f"| {loft} | {len(rows)} ({glabel}) "
              f"| {L.min():.1f} / {Lm:.1f} / {L.max():.1f} | {area:.0f} | {w:.1f} / {anchor:.1f} / {tissue:.1f} | {A:.2f} "
              f"| {' / '.join(f'{v:.2f}' for v in Fp)} | {' / '.join(f'{v:.3f}' for v in sp)} "
              f"({' / '.join(f'{p / Lm:.2f}' for p in pts)}) | {t1[0][0]:.4g}, {t1[0][1]:.2f} ({100 * e1:.0f} %) "
              f"| {t2[0][0]:.3g}, {t2[0][1]:.2f} + {t2[1][0]:.3g}, {t2[1][1]:.2f} ({100 * e2:.0f} %) "
              f"| {s['c']:.4g}, {s['m']:.2f} |")
        # the fits against the section curve at a few strains
        for p in (10.0, 25.0, 36.5, 43.5, umax):
            ep = p / Lm
            lam = 1 + ep
            sec = float(np.interp(p, u, sig))
            o1 = sum(ogden_P(lam, c, m) for c, m in t1)
            o2 = sum(ogden_P(lam, c, m) for c, m in t2)
            os_ = ogden_P(lam, s['c'], s['m'])
            print(f'      u {p:5.1f} mm, strain {ep:.2f}: section {sec:.4f} MPa | 1-term {o1:.4f} | 2-term {o2:.4f} | '
                  f'strip fit {os_:.4f} ({os_ / max(sec, 1e-12):.2f}x)')
        out = os.path.join(DIAG_DIR, f'cl_usl_section_{loft}.txt')
        with open(out, 'w') as fh:
            fh.write(f'** {loft}: sum of its {len(rows)} Abaqus connectors / (average width {w:.3f} mm x t {a.t}); '
                     f'strain = u / {Lm:.2f} mm\n*Uniaxial Test Data\n')
            for ee, ss in zip(eps[::10], sig[::10]):
                fh.write(f'{ss:.6g}, {ee:.6g}\n')


if __name__ == '__main__':
    main()
