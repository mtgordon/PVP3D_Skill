"""What holds the lower canal (2026-10-21; the user: the curves' extra length goes into an arch instead of out through the
prolapse in L379 / L381, while L352's length went out). For RUN at times T:
  - introitus: the mean displacement of the BC-AVW / BC-PVW introitus-edge node sets along OUT (the unit vector from the
    cervix centroid to the introitus centroid at rest: the canal's outward axis), and u_y / u_z;
  - strap (PM_avw_bottom_left/right_sheet_solid): mean axial force ~ sum(element axial stress x volume) / strap length, the
    axis from the strap's curve-edge nodes (shared with the curves) to its far nodes, deformed;
  - AVW-Para <-> PM contact (if present): the total normal force = sum(contact pressure x deformed face area) on its surfaces.
usage: py -3.10 lower_canal_hold.py RUN [RUN ...] [--t 0.6,0.76,0.86,1]"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402
from curve_span import dom_nodes, RUNS  # noqa: E402

STRAP = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
SETS = ('BC-AVW-introitus-edge', 'BC-AVW-introitus-edge-Left', 'BC-PVW-introitus-edge', 'BC-PVW-introitus-edge-Left')


def tet_vol(P, e):
    if len(e) == 6:   # penta6 -> 3 tets
        ts = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
    else:
        ts = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]
    return sum(abs(np.dot(np.cross(P[e[b]] - P[e[a]], P[e[c]] - P[e[a]]), P[e[d]] - P[e[a]])) / 6 for a, b, c, d in ts)


def main(runs, ts):
    for r in runs:
        feb = open(os.path.join(RUNS, r, r + '.feb'), encoding='utf-8').read()
        x = Xplt(os.path.join(RUNS, r, r + '.xplt'))
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        X0 = x.X
        intro = np.array(sorted({idx[int(v)] for s in SETS for m in [re.search(rf'<NodeSet name="{s}">([^<]*)</NodeSet>', feb)] if m
                                 for v in re.findall(r'\d+', m.group(1))}))
        cx = np.array([idx[v] for v in dom_nodes(r, lambda n: n == '_PickedSet346')])
        out_ax = X0[intro].mean(0) - X0[cx].mean(0); out_ax /= np.linalg.norm(out_ax)
        names = [d.get('name') for d in ET.fromstring(feb).find('MeshDomains')]
        curve = {idx[v] for v in dom_nodes(r, lambda n: n.startswith('canal_wrap_hex'))}
        straps = []
        for s in STRAP:
            if s in names:
                k = names.index(s)
                conn = np.array([[idx[v] for v in c] for c in x.domains[k]['conn']])
                straps.append((k, conn))
        surfs = [s for s in ('AVW_Para_solid_PM_primary',) if s in [q.get('name') for q in x.surfaces]]   # primary only: the pair's force once
        conv = [k for k, st in enumerate(x.states) if st[1] == 0]
        tt = np.array([x.states[k][0] for k in conv])
        print(f'== {r[:4]} (last t {tt[-1]:.3f}); outward axis (rest) {np.round(out_ax, 2)}; introitus {len(intro)} nodes; strap domains {len(straps)}; PM contact surfaces {surfs or "none"}')
        for T in ts:
            if T > tt[-1] + 1e-9:
                continue
            k = conv[int(np.argmin(abs(tt - T)))]
            U = x.var(k, 'displacement'); P = X0 + U
            ui = U[intro]
            s = f'  t {x.states[k][0]:.3f}: introitus out {(ui @ out_ax).mean():+6.1f} mm (u_y {ui[:, 1].mean():+5.1f}, u_z {ui[:, 2].mean():+5.1f})'
            if straps:
                S = x.var(k, 'stress')
                F = []
                for dk, conn in straps:
                    sig = S.get(dk + 1, S.get(dk))
                    nodes = np.unique(conn)
                    ce = np.array([n for n in nodes if n in curve]); fa = np.array([n for n in nodes if n not in curve])
                    ax = P[fa].mean(0) - P[ce].mean(0); Ls = np.linalg.norm(ax); ax /= Ls
                    tot = 0.0
                    for e, sv in zip(conn, sig):
                        xx, yy, zz, xy, yz, xz = sv[:6]
                        M = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
                        tot += (ax @ M @ ax) * tet_vol(P, e)
                    F.append(tot / Ls)
                s += f' | strap axial force L / R {F[0]:.2f} / {F[1]:.2f} N' if len(F) == 2 else ''
            if surfs:
                cp = x.var(k, 'contact pressure')
                tot = 0.0
                for sf in x.surfaces:
                    if sf.get('name') not in surfs:
                        continue
                    p = cp.get(sf['id'])
                    if p is None:
                        continue
                    for (fid, q), pv in zip(sf['faces'], p):   # face nodes are 0-based plot indices
                        A = 0.5 * np.linalg.norm(np.cross(P[q[1]] - P[q[0]], P[q[2]] - P[q[0]]))
                        if len(q) == 4:
                            A += 0.5 * np.linalg.norm(np.cross(P[q[2]] - P[q[0]], P[q[3]] - P[q[0]]))
                        tot += float(np.ravel(pv)[0]) * A
                s += f' | AVW-Para <-> PM contact force {tot:.2f} N'
            print(s, flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    ts = [float(v) for v in a[a.index('--t') + 1].split(',')] if '--t' in a else [0.6, 0.76, 0.86, 1.0]
    main([v for i, v in enumerate(a) if not v.startswith('--') and (i == 0 or a[i - 1] != '--t')], ts)
