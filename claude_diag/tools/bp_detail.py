"""Ap / Bp detail for one run (pop_q_bump.py's rules): the Ap node, the Bp node, their rest and deformed coordinates (mm,
model frame) and signed distance to the hymenal plane (+ below), plus the five most dependent upper-PVW nodes and the most
dependent node of the WHOLE PVW lumen surface (for comparison; not Bp).
usage: py -3.10 bp_detail.py RUN [--t 1]"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
args = sys.argv[1:]
tw = None
if '--t' in args:
    i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]
run = args[0]


def nodes_of(root):
    mesh = root.find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    return mesh, X


_, XR = nodes_of(ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot())
P = np.array([XR[i] for i in range(23817, 24213)])
c0 = P.mean(axis=0)
nrm = np.linalg.svd(P - c0)[2][-1]
mesh, X = nodes_of(ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot())
dom = {}
for e in mesh.findall('Elements'):
    for el in e:
        for v in el.text.split(','):
            dom.setdefault(int(v), set()).add(e.get('name').replace('_dist', ''))   # 2026-10-06: a local zone domain <name>_dist counts as <name>
surf = {s.get('name'): {int(v) for f in s for v in f.text.split(',')} for s in mesh.findall('Surface')}
pvw = sorted(n for n in surf['SlidingElastic1Primary'] if '_PickedSet64' in dom.get(n, ()))
cx = [n for n in X if '_PickedSet346' in dom.get(n, ())]
n = nrm if np.dot(np.mean([X[v] for v in cx], axis=0) - c0, nrm) < 0 else -nrm
val = lambda xyz: (xyz - c0) @ n
ap = min([v for v in pvw if abs(X[v][0]) < 2.0], key=lambda v: abs(val(X[v]) + 30))
upper = [v for v in pvw if val(X[v]) <= val(X[ap]) + 1e-9]
x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
conv = converged_state_indices(x, os.path.join(RUNS, run, run + '.log')) or list(range(len(x.states)))
k = conv[-1] if tw is None else min(conv, key=lambda j: abs(x.states[j][0] - tw))
u = x.var(k, 'displacement')
idx = {int(v): i for i, v in enumerate(x.node_ids)}
pos = lambda v: X[v] + u[idx[v]]
f = lambda a: '(' + ', '.join(f'{c:7.2f}' for c in a) + ')'
print(f'{run} at t {x.states[k][0]:.3f}; hymenal plane: point {f(c0)}, unit normal (+ = below the hymen) {f(n)}')
print(f'PVW lumen nodes {len(pvw)}, upper (rest value <= Ap rest {val(X[ap]):+.1f} mm) {len(upper)}')
for lab, v in (('Ap', ap), ('Bp', max(upper, key=lambda v: val(pos(v)))), ('whole-PVW max (not Bp)', max(pvw, key=lambda v: val(pos(v))))):
    print(f'{lab:24s} node {v}: rest {f(X[v])} value {val(X[v]):+6.1f} | deformed {f(pos(v))} value {val(pos(v)):+6.1f}')
print('five most dependent upper-PVW nodes:')
for v in sorted(upper, key=lambda v: -val(pos(v)))[:5]:
    print(f'  node {v}: rest {f(X[v])} value {val(X[v]):+6.1f} | deformed {f(pos(v))} value {val(pos(v)):+6.1f}')
