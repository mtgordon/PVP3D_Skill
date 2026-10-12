"""Self-penetration of the tube models' curves (canal_wrap_hex_avw / _pvw), 2026-10-05.
At each requested time: how many of the curves' outer-surface nodes (build_batch164.outer_faces of the run's own .feb) lie
inside a curve hex that is not theirs or a neighbour's (the 2-ring of hexes around the node is skipped), where they started
(left / right, rest z), and the self-contact surface's pressure / gap (curve_all or curve_outer, whichever the run has).
usage: py -3.10 curve_self_pen.py RUN T1,T2,...   (e.g. curve_self_pen.py L169_wrapfull_endspr_selfall_walls035_P1 0.8,0.9,1)"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline',
                                'scripts'))
sys.path.insert(0, os.path.join(os.path.expanduser('~'), '.claude', 'skills', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt, converged_state_indices  # noqa: E402
from build_batch164 import outer_faces, elems  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
TET = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]   # 5-tet split of a hex8


def main(run, ts):
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    H = np.array(elems(t, 'canal_wrap_hex_avw') + elems(t, 'canal_wrap_hex_pvw'))
    faces, _ = outer_faces(t)
    surf_nodes = sorted({n for f in faces for n in f})
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(RUNS, run, run + '.log'))
    idx = {int(v): k for k, v in enumerate(x.node_ids)}
    Hi = np.vectorize(idx.get)(H)
    names = [s.get('name') for s in x.surfaces]
    sname = next((s for s in ('curve_all', 'curve_outer') if s in names), None)
    node_hex = {}
    for h, e in enumerate(H):
        for n in e:
            node_hex.setdefault(n, set()).add(h)
    for tt in ts:
        k = min(conv, key=lambda q: abs(x.states[q][0] - tt))
        P = x.X + x.var(k, 'displacement')
        hx = P[Hi]
        lo, hi = hx.min(1), hx.max(1)
        inside = []
        for n in surf_nodes:
            p = P[idx[n]]
            cand = np.where(np.all((p >= lo - 1e-9) & (p <= hi + 1e-9), axis=1))[0]
            near = set().union(*[set().union(*[node_hex[m] for m in H[h]]) for h in node_hex[n]])
            for h in cand:
                if h in near:
                    continue
                hit = False
                for a, b, c, d in TET:
                    A, B, C, D = hx[h][[a, b, c, d]]
                    try:
                        lam = np.linalg.solve(np.column_stack([B - A, C - A, D - A]), p - A)
                    except np.linalg.LinAlgError:
                        continue
                    if lam.min() > 1e-6 and lam.sum() < 1 - 1e-6:
                        hit = True
                        break
                if hit:
                    inside.append(n)
                    break
        msg = ''
        if inside:
            Q = np.array([x.X[idx[n]] for n in inside])
            msg = 'left %d / right %d; rest z %.1f..%.1f' % ((Q[:, 0] < 1.4).sum(), (Q[:, 0] >= 1.4).sum(), Q[:, 2].min(),
                                                              Q[:, 2].max())
        cs = ''
        if sname:
            si = names.index(sname) + 1
            pr, gp = x.var(k, 'contact pressure')[si], x.var(k, 'contact gap')[si]
            cs = ' | %s faces with pressure %d, max %.2e MPa' % (sname, (pr > 1e-6).sum(), pr.max())
        print('t %.3f: %d of %d outer-surface nodes inside another part of the curve %s%s'
              % (x.states[k][0], len(inside), len(surf_nodes), msg, cs))


if __name__ == '__main__':
    main(sys.argv[1], [float(v) for v in sys.argv[2].split(',')])
