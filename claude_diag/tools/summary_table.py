"""One markdown row per run for the 2026-09-27 reports: t reached, failed attempts, iterations, wall time (log; a run that
did not end prints its last converged time and "running/stopped"), the LA's ballooning (la_stretch.py's measures: total
area ratio, element area stretch median / p90, edge stretch median) and the POP-Q points Ba / Bp / C as Bump et al. 1996
define them (pop_q_bump.py; before 2026-09-28 Ba / Bp were the tracked 3 cm points, i.e. Aa / Ap, equal to them in every
run checked) (signed distance to the hymenal plane PM_Plane, + below the hymen), all at the last converged state, or at --t.
usage: py -3.10 summary_table.py RUN [RUN ...] [--t 1.0]"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
args = sys.argv[1:]
tw = None
if '--t' in args:
    i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]


def area(P):
    if len(P) == 4:
        return 0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[1]))
    return 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))


def nodes_of(root):
    mesh = root.find('Mesh')
    return mesh, {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}


ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot()
_, XR = nodes_of(ref)
P = np.array([XR[i] for i in range(23817, 24213)])
c0 = P.mean(axis=0)
nrm = np.linalg.svd(P - c0)[2][-1]

print('| run | t | failed | iterations | wall | LA area | element stretch median / p90 | edge stretch | Ba [mm] | Bp [mm] | C [mm] |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for run in args:
    txt = open(os.path.join(RUNS, run, run + '.log'), encoding='latin-1').read()
    g = lambda pat: (re.findall(pat, txt) or [''])[-1]
    fails = len(re.findall(r'------- failed to converge', txt))
    wall = g(r'Total elapsed time \.+ : (\d+:\d+:\d+)')
    its = g(r'Total number of equilibrium iterations \.+ : (\d+)')
    if not its:
        its = str(sum(int(v) for v in re.findall(r'^\s+number of iterations\s+:\s+(\d+)', txt, re.M)))
    end = 'N O R M A L' in txt
    root = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot()
    mesh, X = nodes_of(root)
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    T = np.array([s[0] for s in x.states])
    k = int(np.argmin(np.abs(T - tw))) if tw is not None else len(T) - 1
    U = x.var(k, 'displacement')
    # LA stretch
    els = [[idx[int(v)] for v in e.text.split(',')] for b in mesh.findall('Elements') if b.get('name', '').startswith('LA_')
           for e in b]
    A0 = np.array([area(x.X[c]) for c in els])
    X1 = x.X + U
    r = np.array([area(X1[c]) for c in els]) / A0
    L0 = np.array([np.linalg.norm(x.X[c[i]] - x.X[c[(i + 1) % len(c)]]) for c in els for i in range(len(c))])
    st = np.array([np.linalg.norm(X1[c[i]] - X1[c[(i + 1) % len(c)]]) for c in els for i in range(len(c))]) / L0
    # POP-Q
    dom = {}
    for e in mesh.findall('Elements'):
        for el in e:
            for v in el.text.split(','):
                dom.setdefault(int(v), set()).add(e.get('name'))
    surf = {s.get('name'): {int(v) for f in s for v in f.text.split(',')} for s in mesh.findall('Surface')}
    avw = [n for n in surf['SlidingElastic1Secondary'] if '_PickedSet347' in dom.get(n, ())]
    pvw = [n for n in surf['SlidingElastic1Primary'] if '_PickedSet64' in dom.get(n, ())]
    cx = [n for n in X if '_PickedSet346' in dom.get(n, ())]
    nn = nrm if np.dot(np.mean([X[v] for v in cx], axis=0) - c0, nrm) < 0 else -nrm
    val = lambda xyz: (xyz - c0) @ nn
    mid = lambda ids: [v for v in ids if abs(X[v][0]) < 2.0]
    am, pm = mid(avw), mid(pvw)
    # Bump et al. 1996 (2026-09-28): Aa / Ap = the midline points 3 cm proximal to the hymen at rest, tracked; Ba / Bp =
    # the most dependent point of the upper wall (fornix down to Aa / Ap) at this state, found anew (tools/pop_q_bump.py)
    aa = min(am, key=lambda v: abs(val(X[v]) + 30)); ap = min(pm, key=lambda v: abs(val(X[v]) + 30))
    pos = lambda v: X[v] + U[idx[v]]
    ba = max((v for v in avw if val(X[v]) <= val(X[aa]) + 1e-9), key=lambda v: val(pos(v)))
    bp = max((v for v in pvw if val(X[v]) <= val(X[ap]) + 1e-9), key=lambda v: val(pos(v)))
    C = max(val(pos(v)) for v in cx)
    state = f'{T[k]:.3f}' + ('' if end else ' (not ended)')
    print(f'| `{run}` | {state} | {fails} | {its} | {wall or "-"} | x{(r * A0).sum() / A0.sum():.2f} | '
          f'{np.median(r):.2f} / {np.percentile(r, 90):.2f} | {np.median(st):.2f} | {val(pos(ba)):+.1f} | {val(pos(bp)):+.1f} | '
          f'{C:+.1f} |')
