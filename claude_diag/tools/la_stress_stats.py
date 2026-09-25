"""Von Mises stress distribution per domain group (LA, AVW, PVW) at chosen times, to compare with published stress
plots of the Abaqus model (e.g. the rectocele figure: 0 / 70 / 140 cm H2O = t 0 / 0.5 / 1 on the smooth step).

usage: py -3.10 la_stress_stats.py RUN [RUN ...]
FEBio shell stress in the plot file is the element average through the thickness; an Abaqus 'SNEG' contour is the
bottom surface (plus bending), so compare ranges, not single values.
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from xplt import Xplt
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

RUNS = RUNS_DIR
GROUPS = {'LA': lambda n: n.startswith('LA_'), 'AVW': lambda n: n == '_PickedSet347', 'PVW': lambda n: n == '_PickedSet64'}
TIMES = (0.5, 0.75, 1.0, 1.1)


def vm(s):
    xx, yy, zz, xy, yz, xz = (s[:, i] for i in range(6))
    return np.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2) + 3 * (xy ** 2 + yz ** 2 + xz ** 2))


def dname(x, rid, names):
    d = x.domains[rid - 1] if 0 < rid <= len(x.domains) else {}
    return d.get('name') or (names[rid - 1] if 0 < rid <= len(names) else '')


for run in sys.argv[1:]:
    names = [d.get('name') for d in ET.parse(os.path.join(RUNS, run, run + '.feb')).getroot().find('MeshDomains')]
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    conv = [i for i, s in enumerate(x.states) if s[1] == 0]
    times = np.array([x.states[s][0] for s in conv])
    print(f'== {run} (last converged t = {times[-1]:.4f}); von Mises [MPa] median / p90 / max')
    picks = sorted({int(np.argmin(np.abs(times - t))) for t in TIMES if t <= times[-1] + 0.02} | {len(times) - 1})
    for k in picks:
        st = x.var(conv[k], 'stress')
        cells = []
        for g, sel in GROUPS.items():
            # plot-file region ids are 1-based; use the domain names stored in the plot file
            v = [vm(a) for rid, a in st.items() if sel(dname(x, rid, names))]
            if not v:
                raise SystemExit(f'no domains matched group {g}')
            v = np.concatenate(v)
            cells.append(f'{g} {np.median(v):.3f} / {np.percentile(v, 90):.3f} / {v.max():.3f}')
        load = 3 * times[k] ** 2 - 2 * times[k] ** 3 if times[k] <= 1 else 1.0
        print(f'  t = {times[k]:.3f} ({140 * load:5.1f} cm H2O) | ' + ' | '.join(cells))
