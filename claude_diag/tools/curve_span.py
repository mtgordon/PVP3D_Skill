"""Distal-curve span and perineal-body rotation (2026-10-20 night; the hand-off asked for a small tool for the README row "PeB
rotation vs the curve collapse"). For RUN at times T, per side (rest x < 0 = L, > 0 = R): the distance between the centroid of the
curve nodes at rest z ZTOP +- 1 mm and the centroid of the curve's end row (rest z < ZEND, on the perineal body); its ratio to
rest; and the PeB (_PickedSet66) best-fit rotation about -x (deg, the wanted sense positive; peb_rotation.py's kabsch).
In L352 the span shortens ~60 % (26.5 -> 11.1 mm, right, t 0 -> 0.8) while the rotation stalls at ~28 deg from t 0.6.
usage: py -3.10 curve_span.py RUN [RUN ...] [--t 0.5,0.6,0.7,0.8] [--ztop -15] [--zend -39.5]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402
from peb_rotation import kabsch  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def dom_nodes(run, pick):
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    out = set()
    for m in re.finditer(r'\t\t<Elements type="\w+" name="([^"]+)">\n(.*?)\t\t</Elements>\n', t, re.S):
        if pick(m.group(1)):
            out |= {int(v) for c in re.findall(r'>([0-9,]+)</elem>', m.group(2)) for v in c.split(',')}
    return out


def main(runs, ts, ztop, zend):
    for r in runs:
        x = Xplt(os.path.join(RUNS, r, r + '.xplt'))
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        X0 = x.X
        cn = np.array([idx[v] for v in dom_nodes(r, lambda n: n.startswith('canal_wrap_hex') and 'cervix' not in n)])
        peb = np.array([idx[v] for v in dom_nodes(r, lambda n: n == '_PickedSet66')])
        rows = {}
        for s, sg in (('L', -1), ('R', 1)):
            c = cn[np.sign(X0[cn, 0]) == sg]
            rows[s] = (c[abs(X0[c, 2] - ztop) < 1.0], c[X0[c, 2] < zend])
        span0 = {s: np.linalg.norm(X0[a].mean(0) - X0[b].mean(0)) for s, (a, b) in rows.items()}
        conv = [k for k, st in enumerate(x.states) if st[1] == 0] or list(range(len(x.states)))
        tt = np.array([x.states[k][0] for k in conv])
        print(f'== {r} (last t {tt[-1]:.3f}); rest span L / R {span0["L"]:.1f} / {span0["R"]:.1f} mm '
              f'(rows {len(rows["L"][0])}+{len(rows["L"][1])} / {len(rows["R"][0])}+{len(rows["R"][1])} nodes)')
        for T in ts:
            if T > tt[-1] + 1e-6:
                continue
            k = conv[int(np.argmin(abs(tt - T)))]
            P = X0 + x.var(k, 'displacement')
            sp = {s: np.linalg.norm(P[a].mean(0) - P[b].mean(0)) for s, (a, b) in rows.items()}
            _, rv, _ = kabsch(X0[peb], P[peb])
            print(f'  t {x.states[k][0]:.3f}: span L {sp["L"]:5.1f} ({sp["L"] / span0["L"]:.2f}) R {sp["R"]:5.1f} '
                  f'({sp["R"] / span0["R"]:.2f}) mm; PeB rot_-x {-np.degrees(rv[0]):+5.1f} deg')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: a[a.index(k_) + 1] if k_ in a else d   # noqa: E731
    ts = [float(v) for v in g('--t', '0.5,0.6,0.7,0.8').split(',')]
    opts = {'--t', '--ztop', '--zend'}
    runs = [v for i, v in enumerate(a) if not v.startswith('--') and (i == 0 or a[i - 1] not in opts)]
    main(runs, ts, float(g('--ztop', '-15')), float(g('--zend', '-39.5')))
