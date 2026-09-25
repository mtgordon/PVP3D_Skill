"""Axial strain / force in the truss-spring chains at a run's last converged states (is the chain slack?).

usage: py -3.10 chain_force.py RUN [N_STATES]   (default: Posterior_Arcus_Right_springs, last 3 converged states)
Force from the Posterior_Arcus-Hyper / ATLA-Hyper curve (area 1 mm^2, linear extrapolation below 0, as in FEBio).
Also prints each interior chain node's transverse offset from the straight line through its neighbours.
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from xplt import Xplt
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

RUNS = RUNS_DIR
EPS = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.])
F = np.array([0, 1.38501, 3.18666, 5.5303, 8.57895, 12.5448, 17.7036, 24.4143, 33.1437, 44.4993, 59.271])


def force(e):
    return np.interp(e, EPS, F) if e >= 0 else e * F[1] / EPS[1]


def main():
    run = sys.argv[1]
    ns = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    sets = sys.argv[3:] or ['Posterior_Arcus_Right_springs']
    d = os.path.join(RUNS, run)
    mesh = ET.parse(os.path.join(d, run + '.feb')).getroot().find('Mesh')
    x = Xplt(os.path.join(d, run + '.xplt'))
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    conv = [i for i, s in enumerate(x.states) if s[1] == 0][-ns:]
    for name in sets:
        ds = next(s for s in mesh.findall('DiscreteSet') if s.get('name') == name)
        pairs = [tuple(map(int, e.text.split(','))) for e in ds.findall('delem')]
        print(f'== {name}: spring  L0 [mm] | strain / force [N] at t = '
              + ', '.join('%.4f' % x.states[s][0] for s in conv))
        for a, b in pairs:
            ia, ib = idx[a], idx[b]
            L0 = np.linalg.norm(x.X[ia] - x.X[ib])
            cells = []
            for s in conv:
                u = x.var(s, 'displacement')
                e = np.linalg.norm(x.X[ia] + u[ia] - x.X[ib] - u[ib]) / L0 - 1
                cells.append(f'{e:+.3f}/{force(e):6.2f}')
            print(f'{a:6d}-{b:6d} {L0:5.2f} | ' + '  '.join(cells))
        chain = [pairs[0][0]] + [b for _, b in pairs]
        u = x.var(conv[-1], 'displacement')
        P = np.array([x.X[idx[n]] + u[idx[n]] for n in chain])
        off = []
        for k in range(1, len(chain) - 1):
            v = P[k + 1] - P[k - 1]
            w = P[k] - P[k - 1]
            off.append(np.linalg.norm(w - (w @ v) / (v @ v) * v))
        print('transverse offset of interior nodes at last state [mm]: '
              + ' '.join(f'{n}:{o:.2f}' for n, o in zip(chain[1:-1], off)))


if __name__ == '__main__':
    main()
