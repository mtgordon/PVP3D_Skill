"""The fibre-lofts line (the user's request, 2026-09-28): every family ever lofted (the 14 *_fan shell domains of
L9_fitall_edge) as a weak matrix + fibres along the source connector lines that carry the family's connector law.
Strip model (loft_fit_all.py): connector i (length L_i) owns a strip of width w_i (sum w_i L_i = the loft area) and the
loft thickness t; the loft pull at end displacement u is sum_i P(1 + u/L_i) w_i t, fitted to sum_i F_i(u) (each
connector's own force-first table, tension only, constant past the last point). P is the fibre's nominal stress in
uniaxial tension along it, P = 2 lam dW/dI_n (the same for FEBio's coupled and uncoupled fibres):
  fiber-exp-pow (lam0 1):  P = 2 lam ksi (lam^2-1)^(beta-1) exp(alpha (lam^2-1)^beta)
  fiber-pow-linear:        toe (lam < lam0) P = 2 lam ksi (lam^2-1)^(beta-1), ksi = E/4/(beta-1) I0^-1.5 (I0-1)^(2-beta);
                           linear P = 2 b lam - E, b = ksi (I0-1)^(beta-1) + E/2/lam0      (FEBioMech source, master)
Not fiber-CDF: in FEBio 4.13 its stress is (E/2) lam (lam^2-1) int f, inconsistent with its own tangent
(fibre_mini/fibre_cdf_mini.py).
usage: py -3.10 fibre_lofts.py [--insitu] [--fit]"""
import os
import sys
from collections import OrderedDict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from paths import RUNS_DIR  # noqa: E402
from loft_survey import survey, tri_area  # noqa: E402
from loft_fit_all import strip_widths, table_force  # noqa: E402

L9 = os.path.join(RUNS_DIR, 'L9_fitall_edge', 'L9_fitall_edge.feb')
SPRINGS = 'L87_springs_newline_rhoi0'
T_LOFT = 0.49
# loft -> the springs line's DiscreteSet names (variants6 naming; a trailing '_' = every numbered set; P-arcus: one set
# for both sides)
SETS = {'AVW-Para-L_fan': ['AVW-Para-L_conn'], 'AVW-Para-R_fan': ['AVW-Para-R_conn'],
        'CL-L_fan': ['CL-L_conn_'], 'CL-R_fan': ['CL-R_conn_'], 'USL-L_fan': ['USL-L_conn_'], 'USL-R_fan': ['USL-R_conn_'],
        'PM_fan': ['PM_conn'], 'PM_PeB_Left_fan': ['PM_PeB_Left_conn'], 'PM_PeB_Right_fan': ['PM_PeB_Right_conn'],
        'PM_avw_bottom_left_fan': ['PM_avw_bottom_left_conn'], 'PM_avw_bottom_right_fan': ['PM_avw_bottom_right_conn'],
        'PeB-constrin_fan': ['PeB-constrin_conn'], 'P-arcus-L_fan': ['Parcus_conn'], 'P-arcus-R_fan': ['Parcus_conn']}


def set_names(loft, names):
    return [s for s in names if any(s == p or (p.endswith('_') and s.startswith(p)) for p in SETS[loft])]


def families(feb=L9):
    """{loft: dict(rows, L, w, area, tabs, order)} for the 14 lofts, connectors as loft_fit_all.fit_all selects them."""
    f, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(feb)
    by_loft = OrderedDict()
    for fam, (loft, rows) in res.items():
        if loft:
            by_loft.setdefault(loft, []).extend(rows)
    sprung = {frozenset(p) for pairs in f.discsets.values() for p in pairs}
    out = OrderedDict()
    for loft, rows in by_loft.items():
        lnodes = sorted({n for c in f.elem_blocks[loft][1].values() for n in c})
        X = np.array([f.nodes[n] for n in lnodes])
        on = [r for r in rows if np.linalg.norm(X - r['xb'], axis=1).min() < 0.01
              and frozenset((r['fa'], r['fb'])) not in sprung]
        area = sum(tri_area(f, c) for c in f.elem_blocks[loft][1].values())
        L, w, scale, order = strip_widths(on, area)
        tabs = [np.array(behav[r['behavior']][0]['table']) for r in on]
        out[loft] = dict(rows=on, L=L, w=w, area=area, tabs=tabs, order=order, feb=f,
                         left=[r['elset'] for r in rows if r not in on])
    return out


def P_exppow(lam, ksi, alpha, beta):
    x = np.clip(lam ** 2 - 1, 0, None)
    return 2 * lam * ksi * x ** (beta - 1) * np.exp(alpha * x ** beta)


def P_powlin(lam, E, beta, lam0):
    I0 = lam0 ** 2
    ksi = E / 4 / (beta - 1) * I0 ** -1.5 * (I0 - 1) ** (2 - beta)
    b = ksi * (I0 - 1) ** (beta - 1) + E / 2 / lam0
    x = np.clip(lam ** 2 - 1, 0, None)
    return np.where(lam < lam0, 2 * lam * ksi * x ** (beta - 1), 2 * b * lam - E)


LAWS = {'exp-pow': P_exppow, 'pow-linear': P_powlin}


def target(fam, u):
    return np.stack([table_force(tb, u) for tb in fam['tabs']], axis=1).sum(axis=1)


def pull(fam, law, u, *p):
    lam = 1 + u[:, None] / fam['L'][None, :]
    return (law(lam, *p) * fam['w'][None, :] * T_LOFT).sum(axis=1)


def fit(fam, u):
    """Best fiber-exp-pow and fiber-pow-linear by relative least squares over u (scale closed-form, shape on a grid)."""
    T = target(fam, u)
    ok = T > 1e-3 * T.max()
    res = {}

    def best(law, grid):
        bb = None
        with np.errstate(over='ignore', invalid='ignore'):
            for shape in grid:
                g = pull(fam, law, u, 1.0, *shape)
                if not np.all(np.isfinite(g[ok])) or g[ok].max() <= 0:
                    continue
                c = np.sum(g[ok] / T[ok]) / np.sum(g[ok] ** 2 / T[ok] ** 2)
                err = np.sqrt(np.mean(((c * g[ok] - T[ok]) / T[ok]) ** 2))
                if bb is None or err < bb[0]:
                    bb = (err, (c,) + tuple(shape))
        return bb
    res['exp-pow'] = best(P_exppow, [(a, b) for a in np.r_[0, np.geomspace(0.01, 300, 80)]
                                     for b in (2, 2.25, 2.5, 3, 3.5, 4)])
    res['pow-linear'] = best(P_powlin, [(b, l0) for b in (2, 2.5, 3, 4)
                                        for l0 in np.r_[1.005, 1.01, 1.02, np.arange(1.03, 2.0, 0.02)]])
    return T, res


LIN = dict(beta=2.0, lam0=1.005)       # the linear fibre: fiber-pow-linear, P ~ E (lam - 1) from lam 1.005


def P_lin(lam, E):
    return P_powlin(lam, E, LIN['beta'], LIN['lam0'])


def fit2(fam, u, w=None):
    """Two fibres along the same direction: fiber-pow-linear (E1; beta 2, lam0 1.005, i.e. P ~ E1 (lam - 1), the
    tables' first segment) + fiber-exp-pow (ksi2, alpha, beta; lam0 1; the J-shaped stiffening). (E1, ksi2) >= 0 by
    least squares on the relative error for each (alpha, beta) on a grid. Returns (rms, E1, ksi2, alpha, beta)."""
    T = target(fam, u)
    ok = T > 1e-3 * T.max()
    wt = np.ones_like(u) if w is None else w
    g1 = pull(fam, P_lin, u, 1.0)
    best = None
    with np.errstate(over='ignore', invalid='ignore'):
        for beta in (2.0, 2.5, 3.0, 4.0):
            for alpha in np.r_[0, np.geomspace(0.01, 500, 90)]:
                g2 = pull(fam, P_exppow, u, 1.0, alpha, beta)
                if not np.all(np.isfinite(g2[ok])):
                    continue
                A = np.stack([g1[ok], g2[ok]], axis=1) / T[ok, None] * np.sqrt(wt[ok, None])
                bvec = np.sqrt(wt[ok])
                sol = None
                for cols in ((0, 1), (0,), (1,)):
                    x, *_ = np.linalg.lstsq(A[:, cols], bvec, rcond=None)
                    if np.all(x >= 0):
                        full = np.zeros(2); full[list(cols)] = x
                        r = np.sqrt(np.mean(((A @ full - bvec) / np.sqrt(wt[ok])) ** 2))
                        if sol is None or r < sol[0]:
                            sol = (r, full)
                if sol and (best is None or sol[0] < best[0]):
                    best = (sol[0], sol[1][0], sol[1][1], alpha, beta)
    return T, best


def P_two(lam, E1, ksi2, alpha, beta):
    return P_lin(lam, E1) + P_exppow(lam, ksi2, alpha, beta)


def seg_dist(p, a, b):
    ab = b - a
    s = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0, 1)
    return np.linalg.norm(p - (a + s * ab))


def fibre_dirs(fam, loft):
    """Per element of the loft (file order): the unit fibre direction = the inverse-distance blend of the two nearest
    connector lines' directions (anchor -> tissue end), projected into the element's plane; also the element normal and
    how far out of the plane the blended direction was [deg]. Returns (eids, A, N, out_of_plane_deg)."""
    f = fam['feb']
    et, d = f.elem_blocks[loft]
    XA = np.array([r['xa'] for r in fam['rows']])
    XB = np.array([r['xb'] for r in fam['rows']])
    Tl = XB - XA
    Tl /= np.linalg.norm(Tl, axis=1)[:, None]
    eids, A, N, off = [], [], [], []
    for eid, conn in d.items():
        P = np.array([f.nodes[n] for n in conn[:3]])
        c = P.mean(axis=0)
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n /= np.linalg.norm(n)
        dist = np.array([seg_dist(c, a, b) for a, b in zip(XA, XB)])
        k = np.argsort(dist)[:2]
        v = Tl[k[0]].copy()
        if len(k) > 1:
            t2 = Tl[k[1]] * np.sign(np.dot(Tl[k[1]], Tl[k[0]]) or 1.0)
            w1, w2 = 1 / (dist[k[0]] + 1e-6), 1 / (dist[k[1]] + 1e-6)
            v = (w1 * v + w2 * t2) / (w1 + w2)
        v /= np.linalg.norm(v)
        vp = v - np.dot(v, n) * n
        off.append(np.degrees(np.arcsin(min(1.0, abs(np.dot(v, n))))))
        vp /= np.linalg.norm(vp)
        eids.append(eid); A.append(vp); N.append(n)
    return eids, np.array(A), np.array(N), np.array(off)


def insitu(run=SPRINGS, t=None):
    """(time, {set name: elongations [mm] of its springs}) at the last converged state (or the one nearest t)."""
    from febmodel import Feb
    from xplt_reader import Xplt, converged_state_indices
    d = os.path.join(RUNS_DIR, run)
    f = Feb(os.path.join(d, run + '.feb'))
    x = Xplt(os.path.join(d, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(d, run + '.log')) or list(range(len(x.states)))
    k = conv[-1] if t is None else min(conv, key=lambda j: abs(x.states[j][0] - t))
    u = x.var(k, 'displacement')
    idx = {int(v): i for i, v in enumerate(x.node_ids)}
    out = {}
    for nm, prs in f.discsets.items():
        e = []
        for a, b in prs:
            X = f.nodes[b] - f.nodes[a]
            e.append(np.linalg.norm(X + u[idx[b]] - u[idx[a]]) - np.linalg.norm(X))
        out[nm] = np.array(e)
    return x.states[k][0], out


def pair_elong(run, pairs, t=None):
    """(time, elongations [mm] of the given (a, b) node pairs) in run at its last converged state (or nearest t)."""
    from febmodel import Feb
    from xplt_reader import Xplt, converged_state_indices
    d = os.path.join(RUNS_DIR, run)
    f = Feb(os.path.join(d, run + '.feb'))
    x = Xplt(os.path.join(d, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(d, run + '.log')) or list(range(len(x.states)))
    k = conv[-1] if t is None else min(conv, key=lambda j: abs(x.states[j][0] - t))
    u = x.var(k, 'displacement')
    idx = {int(v): i for i, v in enumerate(x.node_ids)}
    e = [np.linalg.norm(f.nodes[b] - f.nodes[a] + u[idx[b]] - u[idx[a]]) - np.linalg.norm(f.nodes[b] - f.nodes[a])
         for a, b in pairs]
    return x.states[k][0], np.array(e)


def loft_strain(run, lofts, t=None):
    """Per loft element in run: the stretch along its fibre (lam_f), the in-plane principal stretches and the angle of
    the largest from the fibre, against the stretch of its nearest connector line (lam_c = current / rest length of
    the connector's node pair). Prints one row per loft: how the loft deforms compared with the springs it replaced."""
    from febmodel import Feb
    from xplt_reader import Xplt, converged_state_indices
    d = os.path.join(RUNS_DIR, run)
    f = Feb(os.path.join(d, run + '.feb'))
    x = Xplt(os.path.join(d, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(d, run + '.log')) or list(range(len(x.states)))
    k = conv[-1] if t is None else min(conv, key=lambda j: abs(x.states[j][0] - t))
    u = x.var(k, 'displacement')
    idx = {int(v): i for i, v in enumerate(x.node_ids)}
    pos = lambda n: f.nodes[n] + u[idx[n]]
    fams = families()
    print(f'{run} t {x.states[k][0]:.3f}')
    print('| loft | elements | lam_c median (connectors) | lam_f median / p10 / p90 | lam_f / lam_c median | '
          'max principal / lam_f median | principal > 30 deg off the fibre | min principal < 1 |')
    print('|---|---|---|---|---|---|---|---|')
    for loft in lofts:
        fam = fams[loft]
        eids, A, N, off = fibre_dirs(fam, loft)
        rows = fam['rows']
        XA = np.array([r['xa'] for r in rows]); XB = np.array([r['xb'] for r in rows])
        lam_c_rows = np.array([np.linalg.norm(pos(r['fb']) - pos(r['fa'])) / np.linalg.norm(f.nodes[r['fb']] - f.nodes[r['fa']])
                               for r in rows])
        lf, lc, pr, ang, mn = [], [], [], [], []
        conn = f.elem_blocks[loft][1]
        for eid, a in zip(eids, A):
            c = conn[eid]
            X = np.array([f.nodes[n] for n in c[:3]]); xx = np.array([pos(n) for n in c[:3]])
            t1 = a; nrm = np.cross(X[1] - X[0], X[2] - X[0]); nrm /= np.linalg.norm(nrm); t2 = np.cross(nrm, t1)
            R2 = np.array([[np.dot(X[i] - X[0], t1), np.dot(X[i] - X[0], t2)] for i in (1, 2)]).T   # 2x2 reference
            Fm = np.stack([xx[1] - xx[0], xx[2] - xx[0]], axis=1) @ np.linalg.inv(R2)             # 3x2
            C = Fm.T @ Fm
            w, V = np.linalg.eigh(C)
            lam_f = np.sqrt(C[0, 0])
            cen = X.mean(axis=0)
            j = int(np.argmin([seg_dist(cen, p, q) for p, q in zip(XA, XB)]))
            lf.append(lam_f); lc.append(lam_c_rows[j]); pr.append(np.sqrt(w[-1]))
            ang.append(np.degrees(np.arccos(min(1.0, abs(V[0, -1])))))
            mn.append(np.sqrt(max(w[0], 0)))
        lf, lc, pr, ang, mn = map(np.array, (lf, lc, pr, ang, mn))
        print(f'| {loft} | {len(lf)} | {np.median(lam_c_rows):.3f} | {np.median(lf):.3f} / {np.percentile(lf, 10):.3f} / '
              f'{np.percentile(lf, 90):.3f} | {np.median(lf / lc):.2f} | {np.median(pr / lf):.2f} | '
              f'{100 * np.mean(ang > 30):.0f} % | {100 * np.mean(mn < 1):.0f} % |')


def compare(run, ref=SPRINGS, t=None):
    """The in-situ check (SKILL item 12): each family's connector-end elongations (the springs line's node pairs) in run
    against ref at run's last converged t (or t). Prints one markdown row per family."""
    from febmodel import Feb
    fr = Feb(os.path.join(RUNS_DIR, ref, ref + '.feb'))
    tr, _ = pair_elong(run, [], t)
    print(f'| family | pairs | {ref} median / max [mm] | {run} median / max [mm] | ratio of medians |')
    print('|---|---|---|---|---|')
    done = set()
    for loft in SETS:
        sets = set_names(loft, fr.discsets)
        key = tuple(sets)
        if not sets or key in done:
            continue
        done.add(key)
        pairs = [p for s in sets for p in fr.discsets[s]]
        _, e0 = pair_elong(ref, pairs, tr)
        _, e1 = pair_elong(run, pairs, tr)
        r = np.median(e1) / np.median(e0) if abs(np.median(e0)) > 1e-9 else np.nan
        print(f'| {", ".join(sets)} | {len(pairs)} | {np.median(e0):+.2f} / {e0.max():+.2f} | {np.median(e1):+.2f} / '
              f'{e1.max():+.2f} | {r:.2f} |')
    return tr


if __name__ == '__main__':
    if '--compare' in sys.argv:
        i = sys.argv.index('--compare')
        tt = float(sys.argv[i + 2]) if len(sys.argv) > i + 2 else None
        print('t', compare(sys.argv[i + 1], t=tt))
        sys.exit()
    fams = families()
    t, el = insitu() if ('--insitu' in sys.argv or '--fit' in sys.argv) else (None, {})
    for loft, fam in fams.items():
        sets = set_names(loft, el)
        e = np.concatenate([el[s] for s in sets]) if sets else np.array([np.nan])
        tabs = {tuple(map(tuple, tb)) for tb in fam['tabs']}
        print(f"\n== {loft}: {len(fam['rows'])} connectors ({len(tabs)} tables), L {fam['L'].min():.1f}-{fam['L'].max():.1f} mm, "
              f"area {fam['area']:.0f} mm2, strip w {fam['w'].min():.2f}-{fam['w'].max():.2f} (sum {fam['w'].sum():.1f}) mm"
              + (f"; left out {fam['left']}" if fam['left'] else ''))
        if sets:
            print(f"   springs line t {t:.2f}, sets {sets}: elongation min / median / max "
                  f"{e.min():+.2f} / {np.median(e):+.2f} / {e.max():+.2f} mm")
        for tb in tabs:
            print('   table (F, u):', ' '.join(f'({a:.4g},{b:.3g})' for a, b in tb))
        if '--fit' in sys.argv:
            umax = max(2.0, 1.5 * np.nanmax(e)) if sets else 20.0
            u = np.linspace(0.02 * umax, umax, 60)
            T, res = fit(fam, u)
            for law, r in res.items():
                if r is None:
                    print(f'   fit {law}: none'); continue
                err, p = r
                P = pull(fam, LAWS[law], u, *p)
                pts = ' | '.join(f'u {u[j]:.1f}: {T[j]:.3g} vs {P[j]:.3g}' for j in (5, 20, 40, 59))
                print(f"   fit {law:10s} u <= {umax:.1f} mm: {', '.join(f'{v:.4g}' for v in p)}; rms {100 * err:.1f} % "
                      f"| N: {pts}")
