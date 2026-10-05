"""POP-Q points as Bump et al. 1996 define them (Am J Obstet Gynecol 175:10-17; the user's check, 2026-09-28). The plane of
the hymen is the fixed reference: negative = above (proximal to) it, positive = below (distal to) it.
  Aa / Ap  a point fixed in the midline of the anterior / posterior wall 3 cm proximal to the hymen (Bump: Aa 3 cm from the
           external urethral meatus; the model has no urethra, so 3 cm from the hymen as for Ap), tracked; -3..+3 cm
  Ba / Bp  the most distal (most dependent) position of ANY part of the upper anterior / posterior wall, from the fornix
           down to Aa / Ap; -3 cm without prolapse
  C        the most distal edge of the cervix
The walls are the canal's lumen surfaces (AVW: SlidingElastic1 secondary, AVW nodes; PVW: SlidingElastic1 primary, PVW
nodes); "upper" = the nodes whose rest value is at or proximal to Aa's / Ap's. Plane and signs as tools/pop_q.py (the
source's hymenal reference PM_Plane). tools/pop_q.py's "Ba / Bp" are the tracked points, i.e. Aa / Ap here; Luo et al. 2015
describe their Ba / Bp as "the location on the anterior (Ba) and posterior (Bp) vaginal walls that lie 3 cm above the
hymenal ring" while citing Bump. Both are printed.
usage: py -3.10 pop_q_bump.py RUN [RUN ...] [--t 0.85]   (default: the last converged state)
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
args = sys.argv[1:]
tw = None
if '--t' in args:
    i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]


def nodes_of(root):
    mesh = root.find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    return mesh, X


ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot()
_, XR = nodes_of(ref)
P = np.array([XR[i] for i in range(23817, 24213)])
c0 = P.mean(axis=0)
nrm = np.linalg.svd(P - c0)[2][-1]
print('| run | t | Aa | Ba | Ap | Bp | C [mm] | Ba / Bp leading node (rest value, x) |')
print('|---|---|---|---|---|---|---|---|')
for run in args:
    mesh, X = nodes_of(ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot())
    dom = {}
    for e in mesh.findall('Elements'):
        for el in e:
            for v in el.text.split(','):
                dom.setdefault(int(v), set()).add(e.get('name'))
    surf = {s.get('name'): {int(v) for f in s for v in f.text.split(',')} for s in mesh.findall('Surface')}
    avw = sorted(n for n in surf['SlidingElastic1Secondary'] if '_PickedSet347' in dom.get(n, ()))
    pvw = sorted(n for n in surf['SlidingElastic1Primary'] if '_PickedSet64' in dom.get(n, ()))
    cx = [n for n in X if '_PickedSet346' in dom.get(n, ())]
    n = nrm if np.dot(np.mean([X[v] for v in cx], axis=0) - c0, nrm) < 0 else -nrm   # cervix inside = negative
    val = lambda xyz: (xyz - c0) @ n
    mid = lambda ids: [v for v in ids if abs(X[v][0]) < 2.0]
    aa = min(mid(avw), key=lambda v: abs(val(X[v]) + 30))
    ap = min(mid(pvw), key=lambda v: abs(val(X[v]) + 30))
    upper_a = [v for v in avw if val(X[v]) <= val(X[aa]) + 1e-9]
    upper_p = [v for v in pvw if val(X[v]) <= val(X[ap]) + 1e-9]
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    conv = converged_state_indices(x, os.path.join(RUNS, run, run + '.log')) or list(range(len(x.states)))
    k = conv[-1] if tw is None else min(conv, key=lambda j: abs(x.states[j][0] - tw))
    u = x.var(k, 'displacement')
    idx = {int(v): i for i, v in enumerate(x.node_ids)}
    pos = lambda v: X[v] + u[idx[v]]
    ba = max(upper_a, key=lambda v: val(pos(v)))
    bp = max(upper_p, key=lambda v: val(pos(v)))
    lead = f'{ba} ({val(X[ba]):+.0f}, x {X[ba][0]:+.0f}) / {bp} ({val(X[bp]):+.0f}, x {X[bp][0]:+.0f})'
    print(f'| `{run}` | {x.states[k][0]:.3f} | {val(pos(aa)):+.1f} | {val(pos(ba)):+.1f} | {val(pos(ap)):+.1f} | '
          f'{val(pos(bp)):+.1f} | {max(val(pos(v)) for v in cx):+.1f} | {lead} |')
