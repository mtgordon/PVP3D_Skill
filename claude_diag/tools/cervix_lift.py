"""Cervix lift and drop (2026-10-21; with the curve pressure on, L359's cervix rises ~22 mm above rest by t 0.6-0.86, then drops
42 mm at t 0.88-0.93 before the run crashes: "the vaginal wall flies up and then slams back down"). For RUN at times T: the
cervix (_PickedSet346) mean displacement u_y / u_z (mm) and, with --pvw, the PVW (_PickedSet64) hex volume ratio range.
usage: py -3.10 cervix_lift.py RUN [RUN ...] [--t 0.4,0.6,0.8,0.86,0.9,0.92,last]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402
from curve_span import dom_nodes, RUNS  # noqa: E402


def main(runs, ts):
    for r in runs:
        x = Xplt(os.path.join(RUNS, r, r + '.xplt'))
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        cx = np.array([idx[v] for v in dom_nodes(r, lambda n: n == '_PickedSet346')])
        conv = [k for k, s in enumerate(x.states) if s[1] == 0] or list(range(len(x.states)))
        tt = np.array([x.states[k][0] for k in conv])
        out, seen = [], set()
        for T in ts:
            T = tt[-1] if T == 'last' else float(T)
            if T > tt[-1] + 1e-9:
                continue
            k = conv[int(np.argmin(abs(tt - T)))]
            if k in seen:
                continue
            seen.add(k)
            U = x.var(k, 'displacement')[cx]
            out.append(f'{x.states[k][0]:.3f}: y {U[:, 1].mean():+5.1f} z {U[:, 2].mean():+5.1f}')
        print(f'{r[:4]} (last t {tt[-1]:.3f}) cervix mean u: ' + ' | '.join(out))


if __name__ == '__main__':
    a = sys.argv[1:]
    ts = a[a.index('--t') + 1].split(',') if '--t' in a else ['0.4', '0.6', '0.8', '0.86', '0.9', '0.92', 'last']
    main([v for i, v in enumerate(a) if not v.startswith('--') and (i == 0 or a[i - 1] != '--t')], ts)
