"""Range-of-motion report (2026-10-09, the user's goal: "the Bp and Ba can get at least to the points that they got to in the
paper. If the paper has apical location check for that too. Healthy case should stay unprolapsed"). For each run: Ba, Bp, C
(Bump, as tools/pop_q_bump.py) at rest and at the end state (or --t), their motion (end - rest, mm along the hymenal normal),
and the target of the run's case, read from the name's last part (build_batch243 cases):
  healthy Ba, Bp < 0 | P1 Bp >= +4 | P2 Bp >= +9 | F7A Bp >= +10 | F8A Ba >= +10 | F6a0 C motion >= 25 | F6a60 C motion >= 35
(Fig. 6 also gives anterior / posterior motion: F6a0 24 / 41 mm, F6a60 29 / 50 mm, printed as a second check.)
usage: py -3.10 rom_report.py RUN [RUN ...] [--t 0.9]"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
TARGET = {'healthy': ('Ba < 0 and Bp < 0', lambda r: r['Ba'] < 0 and r['Bp'] < 0),
          'P1': ('Bp >= +4', lambda r: r['Bp'] >= 4), 'P2': ('Bp >= +9', lambda r: r['Bp'] >= 9),
          'F7A': ('Bp >= +10', lambda r: r['Bp'] >= 10), 'F8A': ('Ba >= +10', lambda r: r['Ba'] >= 10),
          'F6a0': ('C moves >= 25 (Ba 24, Bp 41)', lambda r: r['dC'] >= 25),
          'F6a60': ('C moves >= 35 (Ba 29, Bp 50)', lambda r: r['dC'] >= 35),
          'MAX': ('Ba >= +10, Bp >= +10, C moves >= 35', lambda r: r['Ba'] >= 10 and r['Bp'] >= 10 and r['dC'] >= 35)}


def nodes_of(root):
    mesh = root.find('Mesh')
    return mesh, {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}


def measure(run, tw=None):
    ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot()
    _, XR = nodes_of(ref)
    P = np.array([XR[i] for i in range(23817, 24213)])
    c0 = P.mean(axis=0); nrm = np.linalg.svd(P - c0)[2][-1]
    mesh, X = nodes_of(ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot())
    dom = {}
    for e in mesh.findall('Elements'):
        for el in e:
            for v in el.text.split(','):
                dom.setdefault(int(v), set()).add(e.get('name').replace('_dist', ''))
    surf = {s.get('name'): {int(v) for f in s for v in f.text.split(',')} for s in mesh.findall('Surface')}
    avw = sorted(n for n in surf['SlidingElastic1Secondary'] if '_PickedSet347' in dom.get(n, ()))
    pvw = sorted(n for n in surf['SlidingElastic1Primary'] if '_PickedSet64' in dom.get(n, ()))
    cx = [n for n in X if '_PickedSet346' in dom.get(n, ())]
    nn = nrm if np.dot(np.mean([X[v] for v in cx], axis=0) - c0, nrm) < 0 else -nrm
    val = lambda xyz: (xyz - c0) @ nn  # noqa: E731
    mid = lambda ids: [v for v in ids if abs(X[v][0]) < 2.0]  # noqa: E731
    aa = min(mid(avw), key=lambda v: abs(val(X[v]) + 30)); ap = min(mid(pvw), key=lambda v: abs(val(X[v]) + 30))
    ua = [v for v in avw if val(X[v]) <= val(X[aa]) + 1e-9]; up = [v for v in pvw if val(X[v]) <= val(X[ap]) + 1e-9]
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(RUNS, run, run + '.log')) or list(range(len(x.states)))
    k = conv[-1] if tw is None else min(conv, key=lambda j: abs(x.states[j][0] - tw))
    u = x.var(k, 'displacement'); idx = {int(v): i for i, v in enumerate(x.node_ids)}
    pos = lambda v: X[v] + u[idx[v]]  # noqa: E731
    r = {'t': x.states[k][0],
         'Ba': max(val(pos(v)) for v in ua), 'Bp': max(val(pos(v)) for v in up), 'C': max(val(pos(v)) for v in cx),
         'Ba0': max(val(X[v]) for v in ua), 'Bp0': max(val(X[v]) for v in up), 'C0': max(val(X[v]) for v in cx)}
    for p in ('Ba', 'Bp', 'C'):
        r['d' + p] = r[p] - r[p + '0']
    return r


if __name__ == '__main__':
    args = sys.argv[1:]; tw = None
    if '--t' in args:
        i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]
    print('| run | t | Ba | Bp | C | motion Ba / Bp / C | target | met |')
    print('|---|---|---|---|---|---|---|---|')
    for run in args:
        r = measure(run, tw)
        case = run.rsplit('_', 1)[-1]
        tg, ok = TARGET.get(case, ('', None))
        met = '' if ok is None else ('**yes**' if ok(r) else 'no')
        print(f"| `{run}` | {r['t']:.3f} | {r['Ba']:+.1f} | {r['Bp']:+.1f} | {r['C']:+.1f} | "
              f"{r['dBa']:+.1f} / {r['dBp']:+.1f} / {r['dC']:+.1f} | {tg} | {met} |")
