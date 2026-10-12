"""Node speed check (2026-10-06, for the curve flapping): at each given t, the max and 95th-percentile node displacement rate
in mm per 0.001 t (over the converged states >= 0.002 t back), and the domains of the 50 fastest nodes. A settling run stays
~0.05 p95 (L169); L189 flapped from t ~0.68 (p95 0.14 -> 0.29 at 0.70 -> 2-4).
usage: py -3.10 node_speed.py RUN t1,t2,...   ('last' = the last state)"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
run, ts = sys.argv[1], sys.argv[2].split(',')
mesh = ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot().find('Mesh')
dom = {}
for e in mesh.findall('Elements'):
    for el in e:
        for v in el.text.split(','):
            dom.setdefault(int(v), set()).add(e.get('name'))
x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
T = np.array([s[0] for s in x.states])
for t in ts:
    i = len(T) - 1 if t == 'last' else int(np.argmin(np.abs(T - float(t))))
    j = i - 1
    while j > 0 and T[i] - T[j] < 0.002:
        j -= 1
    d = np.linalg.norm(np.asarray(x.var(i, 'displacement'), float) - np.asarray(x.var(j, 'displacement'), float), axis=1)
    d *= 0.001 / (T[i] - T[j])
    c = {}
    for k in np.argsort(-d)[:50]:
        key = ','.join(sorted(dom.get(int(x.node_ids[k]), {'-'})))
        c[key] = c.get(key, 0) + 1
    print(f'{run} t {T[i]:.3f}: max {d.max():.2f} p95 {np.percentile(d, 95):.3f} mm/0.001t; fastest 50: '
          f'{dict(sorted(c.items(), key=lambda kv: -kv[1])[:3])}')
