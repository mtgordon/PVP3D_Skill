"""POP-Q-style points. NOTE (2026-09-28): the "Ba / Bp" here are the TRACKED points 3 cm above the hymen, i.e. Bump et al.
1996's Aa / Ap; Bump's Ba / Bp (the most dependent point of the upper wall) are in pop_q_bump.py and summary_table.py.
(Luo et al. 2015, J Biomech: Ba / Bp = the points on the anterior / posterior vaginal wall 3 cm above
the hymenal ring; C = the cervix), measured against the hymenal reference plane (the source's display body PM_Plane).
Value = signed distance to the plane in mm, negative above (inside) the hymen, positive below, as in POP-Q.
Ba / Bp: the midline (|x| < 2 mm) nodes of the canal's inner surfaces (AVW: SlidingElastic1 secondary; PVW: SlidingElastic1
primary, PVW nodes only) whose rest value is nearest -30 mm, tracked; also each wall's most dependent midline value
(max) and the cervix's most dependent value. usage: py -3.10 pop_q.py RUN [RUN ...] [--t 0.25,0.5,0.75,1.0]
The plane is taken from L26_springs_la3_pm (PM_Plane nodes 23817-24212; same coordinates in every model)."""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
args = [a for a in sys.argv[1:] if not a.startswith('--')]
ts = [1.0]
if '--t' in sys.argv:
    ts = [float(v) for v in sys.argv[sys.argv.index('--t') + 1].split(',')]
    args = [a for a in args if a != sys.argv[sys.argv.index('--t') + 1]]


def nodes_of(root):
    mesh = root.find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    return mesh, X


ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot()
_, XR = nodes_of(ref)
P = np.array([XR[i] for i in range(23817, 24213)])
c0 = P.mean(axis=0)
nrm = np.linalg.svd(P - c0)[2][-1]
for run in args:
    root = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot()
    mesh, X = nodes_of(root)
    dom = {}
    for e in mesh.findall('Elements'):
        for el in e:
            for v in el.text.split(','):
                dom.setdefault(int(v), set()).add(e.get('name'))
    surf = {s.get('name'): {int(v) for f in s for v in f.text.split(',')} for s in mesh.findall('Surface')}
    avw = [n for n in surf['SlidingElastic1Secondary'] if '_PickedSet347' in dom.get(n, ())]
    pvw = [n for n in surf['SlidingElastic1Primary'] if '_PickedSet64' in dom.get(n, ())]
    cx = [n for n in X if '_PickedSet346' in dom.get(n, ())]
    n = nrm if np.dot(np.mean([X[v] for v in cx], axis=0) - c0, nrm) < 0 else -nrm   # cervix inside = negative
    val = lambda xyz: (xyz - c0) @ n
    mid = lambda ids: [v for v in ids if abs(X[v][0]) < 2.0]
    am, pm = mid(avw), mid(pvw)
    ba = min(am, key=lambda v: abs(val(X[v]) + 30)); bp = min(pm, key=lambda v: abs(val(X[v]) + 30))
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = {int(k): i for i, k in enumerate(x.node_ids)}
    T = np.array([s[0] for s in x.states])
    print(f'{run}: Ba node {ba} (rest {val(X[ba]):+.1f} mm), Bp node {bp} (rest {val(X[bp]):+.1f}); rest most dependent: '
          f'AVW {max(val(X[v]) for v in am):+.1f}, PVW {max(val(X[v]) for v in pm):+.1f}, cervix {max(val(X[v]) for v in cx):+.1f}')
    for tw in ts:
        k = int(np.argmin(np.abs(T - tw)))
        u = x.var(k, 'displacement')
        pos = lambda v: X[v] + u[idx[v]]
        print(f'  t {T[k]:.3f}: Ba {val(pos(ba)):+6.1f}  Bp {val(pos(bp)):+6.1f}  | most dependent: AVW {max(val(pos(v)) for v in am):+6.1f}  '
              f'PVW {max(val(pos(v)) for v in pm):+6.1f}  cervix (C) {max(val(pos(v)) for v in cx):+6.1f} mm')
