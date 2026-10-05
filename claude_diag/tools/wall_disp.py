"""Vaginal-wall movement per region: |u| median / p90 / max and the mean displacement vector (x, y, z), at a few times.

usage: py -3.10 wall_disp.py RUN [RUN ...] [--t 0.5,0.75,1.0]
Regions (the Abaqus VW-PeB part's element sets): AVW (_PickedSet347, Vagina_AVW), PVW (_PickedSet64, Vagina_PVW), cervix /
apex (_PickedSet346, Vagina_Cervix), perineal body (_PickedSet66), and the LA (all LA_* shell domains) for reference.
Each converged state nearest a requested time is used; a node shared by two regions counts in both.
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from xplt import Xplt  # noqa: E402

REGIONS = (('AVW', ('_PickedSet347',)), ('PVW', ('_PickedSet64',)), ('cervix', ('_PickedSet346',)),
           ('PeB', ('_PickedSet66',)), ('LA', ('LA_ICM', 'LA_ICM_tri', 'LA_PCM', 'LA_PCMPRM')))


def main():
    args = sys.argv[1:]
    times = (0.5, 0.75, 1.0)
    if '--t' in args:
        i = args.index('--t')
        times = tuple(float(v) for v in args[i + 1].split(','))
        args = args[:i] + args[i + 2:]
    for run in args:
        names = [d.get('name') for d in ET.parse(os.path.join(RUNS_DIR, run, run + '.feb')).getroot().find('MeshDomains')]
        x = Xplt(os.path.join(RUNS_DIR, run, run + '.xplt'))
        nodes = {}
        for label, doms in REGIONS:
            s = set()
            for i, d in enumerate(x.domains):
                if i < len(names) and names[i] in doms:
                    s.update(int(n) for c in d['conn'] for n in c)
            nodes[label] = np.array(sorted(s))
        conv = [k for k, st in enumerate(x.states) if st[1] == 0]
        tt = np.array([x.states[k][0] for k in conv])
        print(f'== {run} (last converged t = {tt[-1]:.3f})')
        print(f'{"t":>6} {"region":7s} |u| median   p90   max [mm] | mean u (x, y, z) [mm]     | u_y p10 / p50 / p90 / max [mm] (+y = forward)')
        for t in times:
            if t > tt[-1] + 1e-6:
                continue
            k = conv[int(np.argmin(np.abs(tt - t)))]
            u = x.var(k, 'displacement')
            for label, _ in REGIONS:
                un = u[nodes[label]]
                a = np.linalg.norm(un, axis=1)
                mu = un.mean(axis=0)
                uy = un[:, 1]
                print(f'{x.states[k][0]:6.3f} {label:7s} {np.median(a):10.1f} {np.percentile(a, 90):5.1f} {a.max():5.1f}      '
                      f'| ({mu[0]:+5.1f}, {mu[1]:+5.1f}, {mu[2]:+5.1f})     | {np.percentile(uy, 10):+5.1f} / '
                      f'{np.percentile(uy, 50):+5.1f} / {np.percentile(uy, 90):+5.1f} / {uy.max():+5.1f}')


if __name__ == '__main__':
    main()
