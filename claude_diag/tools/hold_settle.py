"""The equilibrium check for a run with its load held past the ramp (batch 28: L14_p3_hold).

usage: py -3.10 hold_settle.py RUN [--t0 1.0] [--every 0.1]
1. Per converged state from t0 - 0.2 on (every ~--every s, plus the last): LA |u| median / p90 / max, the fastest LA
   node and the fastest node of the whole model (speed = |du| / dt between consecutive converged states, or the plotted
   velocity when the run has it) with its domain.
2. The settle fit of the skill (convergence-debugging.md): y = a + b exp(-(t - t0) / tau) to each LA percentile over
   the hold (t >= t0); a is the equilibrium estimate. tau by a grid search, a and b by linear least squares.
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from xplt import Xplt  # noqa: E402


def settle_fit(t, y, t0):
    best = None
    for tau in np.geomspace(0.01, 20.0, 400):
        A = np.column_stack([np.ones_like(t), np.exp(-(t - t0) / tau)])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ coef
        rms = float(np.sqrt(np.mean(r ** 2)))
        if best is None or rms < best[3]:
            best = (coef[0], coef[1], tau, rms)
    return best


def main():
    args = sys.argv[1:]
    t0 = float(args[args.index('--t0') + 1]) if '--t0' in args else 1.0
    every = float(args[args.index('--every') + 1]) if '--every' in args else 0.1
    run = args[0]
    root = ET.parse(os.path.join(RUNS_DIR, run, run + '.feb')).getroot()
    names = [d.get('name') for d in root.find('MeshDomains')]
    x = Xplt(os.path.join(RUNS_DIR, run, run + '.xplt'))
    dom_of = {}
    la = set()
    for i, d in enumerate(x.domains):
        nm = names[i] if i < len(names) else f'dom{i}'
        for c in d['conn']:
            for n in c:
                dom_of.setdefault(int(n), nm)
                if nm.startswith('LA_'):
                    la.add(int(n))
    la_i = np.array(sorted(la))
    has_v = any(it.get('name') == 'velocity' for it in x.dict.get('nodal', []))
    conv = [i for i, s in enumerate(x.states) if s[1] == 0]
    times = np.array([x.states[s][0] for s in conv])
    print(f'== {run}: {len(conv)} converged states, last t = {times[-1]:.4f}; speeds from '
          f'{"the plotted velocity" if has_v else "du/dt between converged states"}')
    print(f'{"t":>6} | LA |u| med   p90   max [mm] | LA max speed | model max speed [mm/s], domain')
    rows = []
    prev_u, prev_t = None, None
    next_print = t0 - 0.2
    for k, s in enumerate(conv):
        t = times[k]
        if t < t0 - 0.25 and k < len(conv) - 1 and times[k + 1] < t0 - 0.2:
            continue
        u = x.var(s, 'displacement')
        un = np.linalg.norm(u[la_i], axis=1)
        p = (np.median(un), np.percentile(un, 90), un.max())
        if t >= t0 - 1e-9:
            rows.append((t, *p))
        if has_v:
            v = np.linalg.norm(x.var(s, 'velocity'), axis=1)
        elif prev_u is not None and t > prev_t:
            v = np.linalg.norm(u - prev_u, axis=1) / (t - prev_t)
        else:
            v = None
        if v is not None and (t >= next_print - 1e-9 or k == len(conv) - 1):
            j = int(np.argmax(v))
            print(f'{t:6.3f} | {p[0]:12.2f} {p[1]:5.2f} {p[2]:5.2f}   | {v[la_i].max():12.1f} | {v[j]:8.1f}, '
                  f'{dom_of.get(int(j), "?")}')
            next_print = t + every
        prev_u, prev_t = u, t
    if len(rows) >= 4:
        R = np.array(rows)
        print(f'settle fit over t = {R[0, 0]:.3f}-{R[-1, 0]:.3f} ({len(R)} states): y = a + b exp(-(t - {t0:g}) / tau)')
        for j, lab in ((1, 'median'), (2, 'p90'), (3, 'max')):
            a, b, tau, rms = settle_fit(R[:, 0], R[:, j], t0)
            print(f'   {lab:6s}: equilibrium a = {a:6.2f} mm (last {R[-1, j]:6.2f}, at t0 {R[0, j]:6.2f}); b = {b:+.2f}, '
                  f'tau = {tau:.3f} s, rms {rms:.3f} mm')


if __name__ == '__main__':
    main()
