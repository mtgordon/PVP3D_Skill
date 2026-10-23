"""Wall midline outward motion by level (2026-10-22; the README 10-21 'wall midline outward' rows as a tool): mean displacement
of the AVW (_PickedSet347) / PVW (_PickedSet64) nodes with |x| < 4 mm at rest along the canal's outward axis (cervix ->
introitus at rest, lower_canal_hold.py), in rest-z bands upper -> lower. A mid-canal bulge (L352: 28 / 29 / 34 / 37 / 23)
is the prolapse; equal values at every level (L379 / L381: 15-19) = the canal moving as a rigid tube (the curve rails).
usage: py -3.10 midline_bulge.py RUN [RUN ...] [--t 0.76,0.95]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402
from curve_span import dom_nodes, RUNS  # noqa: E402
from lower_canal_hold import SETS  # noqa: E402

BANDS = ((0, 15), (-10, 0), (-20, -10), (-30, -20), (-40, -30))


def main(runs, ts):
    for r in runs:
        feb = open(os.path.join(RUNS, r, r + '.feb'), encoding='utf-8').read()
        x = Xplt(os.path.join(RUNS, r, r + '.xplt'))
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        X0 = x.X
        intro = np.array(sorted({idx[int(v)] for s in SETS for m in [re.search(rf'<NodeSet name="{s}">([^<]*)</NodeSet>', feb)] if m
                                 for v in re.findall(r'\d+', m.group(1))}))
        cx = np.array([idx[v] for v in dom_nodes(r, lambda n: n == '_PickedSet346')])
        ax = X0[intro].mean(0) - X0[cx].mean(0); ax /= np.linalg.norm(ax)
        walls = {}
        for w, d in (('AVW', '_PickedSet347'), ('PVW', '_PickedSet64')):
            nd = np.array([idx[v] for v in dom_nodes(r, lambda n, d=d: n == d)])
            nd = nd[abs(X0[nd, 0]) < 4]
            walls[w] = [nd[(X0[nd, 2] >= lo) & (X0[nd, 2] < hi)] for lo, hi in BANDS]
        conv = [k for k, st in enumerate(x.states) if st[1] == 0]
        tt = np.array([x.states[k][0] for k in conv])
        print(f'== {r[:5]} (last t {tt[-1]:.3f}); bands rest z {BANDS} (upper -> lower)')
        for T in ts:
            if T > tt[-1] + 1e-9:
                continue
            k = conv[int(np.argmin(abs(tt - T)))]
            U = x.var(k, 'displacement')
            s = f'  t {x.states[k][0]:.3f}: introitus {(U[intro] @ ax).mean():+5.1f}'
            for w, bs in walls.items():
                s += f' | {w} ' + ' / '.join(f'{(U[b] @ ax).mean():+5.1f}' if len(b) else '  -  ' for b in bs)
            print(s)


if __name__ == '__main__':
    a = sys.argv[1:]
    ts = [0.76, 0.95]
    if '--t' in a:
        i = a.index('--t'); ts = [float(v) for v in a[i + 1].split(',')]; a = a[:i] + a[i + 2:]
    main(a, ts)
