"""LA displacement distribution (not just the peak node) + posterior-arcus chain / connector state per time.

usage: py -3.10 la_disp_stats.py RUN [RUN ...]
Per selected converged time: LA nodal |u| median / p90 / p99 / max and the number of LA nodes above 10, 20 and
30 mm; the posterior-arcus chain strain range; and, where the run has them, the Parcus_conn connectors'
elongation range and total force (Abaqus LA-Y-parcus table, tension-only).
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from xplt import Xplt
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

RUNS = RUNS_DIR
CONN_F = np.array([0, .11667, .28, .56, .84, 1.26, 1.75, 2.42667, 3.26667, 4.41, 5.92667])
CONN_U = np.array([0, 4.7, 9.515, 13.985, 18.8, 23.5, 28.2, 33.015, 37.715, 42.415, 47.115])
TIMES = (0.25, 0.29, 0.33, 0.5, 0.65, 0.8, 0.9, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0)


def run_stats(run):
    root = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot()
    mesh = root.find('Mesh')
    names = [d.get('name') for d in root.find('MeshDomains')]
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    la = set()
    for i, d in enumerate(x.domains):
        if i < len(names) and names[i].startswith('LA_'):
            la.update(int(n) for c in d['conn'] for n in c)
    la_i = np.array(sorted(la))
    sets = {ds.get('name'): [tuple(idx[int(v)] for v in e.text.split(',')) for e in ds.findall('delem')]
            for ds in mesh.findall('DiscreteSet')}
    chain = sets['Posterior_Arcus_Left_springs'] + sets['Posterior_Arcus_Right_springs']
    conn = sets.get('Parcus_conn', [])
    conv = [i for i, s in enumerate(x.states) if s[1] == 0]
    times = np.array([x.states[s][0] for s in conv])
    print(f'== {run}: last converged t = {times[-1]:.4f}')
    hdr = f'{"t":>6} | LA |u| med  p90  p99  max [mm] | n>10 n>20 n>30 | chain strain min/max'
    print(hdr + (' | conn stretch min/max [mm], sum F [N]' if conn else ''))
    picks = sorted({int(np.argmin(np.abs(times - t))) for t in TIMES if t <= times[-1] + 1e-9} | {len(times) - 1})
    for k in picks:
        u = x.var(conv[k], 'displacement')
        un = np.linalg.norm(u[la_i], axis=1)
        p = np.percentile(un, [50, 90, 99])
        P = x.X + u

        def stretch(pairs):
            return np.array([np.linalg.norm(P[a] - P[b]) - np.linalg.norm(x.X[a] - x.X[b]) for a, b in pairs])

        L0 = np.array([np.linalg.norm(x.X[a] - x.X[b]) for a, b in chain])
        e = stretch(chain) / L0
        line = (f'{times[k]:6.3f} | {p[0]:10.1f} {p[1]:4.1f} {p[2]:4.1f} {un.max():4.1f}      | '
                f'{(un > 10).sum():4d} {(un > 20).sum():4d} {(un > 30).sum():4d} | {e.min():+.3f}/{e.max():+.3f}')
        if conn:
            s = stretch(conn)
            F = np.where(s > 0, np.interp(s, CONN_U, CONN_F), 0.0)
            line += f'         | {s.min():+.1f}/{s.max():+.1f}, {F.sum():.1f}'
        print(line)
    print(f'   ({len(la_i)} LA nodes)')


if __name__ == '__main__':
    for r in sys.argv[1:]:
        run_stats(r)
