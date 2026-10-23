"""Batch 332 (2026-10-19; the user, on L328 la60's curve-tip corner node 27586 going in while its neighbours go out: "how about
getting rid of the springs and just use the loft so it's a lot stiffer?"). The curve-tip -> PM fan (PM_avw_bottom L / R sheets,
built by build_batch262 --fan) as a LOFT ONLY: its 5 spring chains a side (curve_end_PM_conn x3 + PM_avw_bottom_conn x2) are
removed and the sheet gets its own 1-term Ogden fitted to those chains (strip model as loft_fit_all.py: connector i = a strip of
its rest length L_i and width w_i = sheet area / sum L, thickness T_SHEET; uniform pull u -> stretch 1 + u / L_i; target = sum of
the chains' own tables, each chain n links in series: F(u) = scale x table(u / n)); c by least squares on relative error over
u in [U0, U1] mm, m on a grid; k = 250 mu0 (skill item 12). The chains' interior nodes stay as sheet nodes.
In la60 the fan chains stretch <= 1.4 mm by t 1 (only the tip-end chain to 27588 / 25900; the rest <= 0.11 mm).
NOT IN SOURCE (a loft in place of connectors; the loft carries compression, the springs did not).
usage: py -3.10 build_batch332.py SRC_RUN NEW_RUN [--u0 0.25] [--u1 4]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch111 import Model111, T_SHEET  # noqa: E402
from variants2 import renumber_discrete_materials  # noqa: E402
from variants import emit  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
M_GRID = np.r_[np.arange(1.5, 10.01, 0.25), np.arange(10.5, 40.01, 0.5)]


def ogden_P(lam, c, m):
    return (c / m) * (lam ** (m - 1) - lam ** (-m / 2 - 1))


def chains(pairs):
    nxt = {a: b for a, b in pairs}
    out = []
    for h in set(nxt) - set(nxt.values()):
        c = [h]
        while c[-1] in nxt:
            c.append(nxt[c[-1]])
        out.append(c)
    return out


def fit(con, area, u0, u1):
    u = np.linspace(u0, u1, 40)
    L = np.array([c[0] for c in con]); w = area / L.sum()
    target = sum(sc * np.interp(u / n, pts[:, 0], pts[:, 1]) for _, n, pts, sc in con)
    best = None
    for m in M_GRID:
        g = sum(ogden_P(1 + u / Li, 1.0, m) * w * T_SHEET for Li, *_ in con)
        c = np.sum(g / target) / np.sum((g / target) ** 2)
        err = np.sqrt(np.mean((c * g / target - 1) ** 2))
        if best is None or err < best[0]:
            best = (err, c, m)
    return best


def main(src, new, u0=0.25, u1=4.0):
    m = Model111(os.path.join(RUNS, src, src + '.feb'))
    X = m.nodes()
    disc = m.root.find('Discrete')
    m.log.append(f'{src} -> the curve-tip fan as a LOFT ONLY (build_batch332.py)')
    for side in ('left', 'right'):
        sheet = f'PM_avw_bottom_{side}_sheet'
        con = []
        for s in (f'curve_end_PM_conn_PM_avw_bottom_{side}_chain', f'PM_avw_bottom_{side}_conn_PM_avw_bottom_{side}_chain'):
            pairs, mat, b, ds = m._spring_set(s)
            pts = np.array([[float(v) for v in p.text.split(',')] for p in mat.find('force').find('points')])
            sc = float(mat.findtext('scale'))
            for c in chains(pairs):
                con.append((sum(np.linalg.norm(X[c[k]] - X[c[k + 1]]) for k in range(len(c) - 1)), len(c) - 1, pts, sc))
            disc.remove(b); disc.remove(mat); m.mesh.remove(ds)
            m.log.append(f'{s}: its {len(pairs)} chain springs removed')
        tris = next(e for e in m.mesh.findall('Elements') if e.get('name') == sheet)
        area = sum(0.5 * np.linalg.norm(np.cross(X[q[1]] - X[q[0]], X[q[2]] - X[q[0]]))
                   for q in ([int(v) for v in e.text.split(',')] for e in tris))
        err, c, mm = fit(con, area, u0, u1)
        mu0 = c / 2
        mat = next(x for x in m.root.find('Material') if x.get('name') == f'{sheet}_mat')
        old = mat.findtext('E')
        for ch in list(mat):
            mat.remove(ch)
        mat.set('type', 'Ogden')
        vals = [('density', '1.06e-09'), ('k', f'{250 * mu0:.6g}'), ('c1', f'{c:.6g}')] + [(f'c{i}', '0') for i in range(2, 7)] \
            + [('m1', f'{mm:g}')] + [(f'm{i}', '1') for i in range(2, 7)]
        for k, v in vals:
            sub = mat.makeelement(k, {})
            sub.text = v
            mat.append(sub)
        m.log.append(f'NOT IN SOURCE: {sheet}_mat neo-Hookean E {old} -> 1-term Ogden c1 {c:.4g}, m1 {mm:g}, k {250 * mu0:.4g} '
                      f'(mu0 {mu0:.4g} MPa, E0 ~{3 * mu0:.4g}), fitted to its {len(con)} removed chains (rest {min(x[0] for x in con):.1f}-'
                      f'{max(x[0] for x in con):.1f} mm, area {area:.0f} mm2, thickness {T_SHEET}) over u {u0}-{u1} mm: rms relative error '
                      f'{err * 100:.1f} %')
    # 2026-10-19 (L333 from L298): spring-patch anchors (build_batch288) left with no spring are removed with their tie rows
    anc = next((b for b in m.mesh.findall('Nodes') if b.get('name') == 'spring_patch_anchors'), None)
    if anc is not None:
        used = {int(v) for ds in m.mesh.findall('DiscreteSet') for e in ds.findall('delem') for v in e.text.split(',')}
        used |= {int(v) for blk in m.mesh.findall('Elements') for e in blk for v in e.text.split(',')}   # L334: loft on anchors
        orphan = {int(n.get('id')) for n in anc if int(n.get('id')) not in used}
        rows = 0
        for c in m.root.iter('constraint'):
            if c.get('name') == 'spring_patch_ties':
                for lc in list(c.findall('linear_constraint')):
                    if int(lc.find('node').get('id')) in orphan:
                        c.remove(lc); rows += 1
        for n in list(anc):
            if int(n.get('id')) in orphan:
                anc.remove(n)
        assert rows == 3 * len(orphan), (rows, len(orphan))
        m.log.append(f'{len(orphan)} spring-patch anchor nodes left without a spring removed ({sorted(orphan)}) with their {rows} '
                     f'spring_patch_ties rows')
    renumber_discrete_materials(m)
    assert not os.path.exists(os.path.join(RUNS, new)), new
    emit(new, m)
    print('\n'.join(m.log))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--u0') + 1]) if '--u0' in a else 0.25, float(a[a.index('--u1') + 1]) if '--u1' in a else 4.0)
