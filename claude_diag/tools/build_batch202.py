"""Batch 202 (2026-10-07, the user away ~4 h: "I wonder if closing up the apical side of the tube would help with that
inversion ... It should be a completely closed tube at the apical side ... the holes are on the corner where the current
curve doesn't go all of the way to the back"). One change on L192_botfull_distpress_cdamp1000_P1 (B: equal PM-arch springs,
distal-half curve pressure, curve damping C 1000; the best tube run), P1:
  L202_distpress_apexclosed_P1   each curve extruded beyond its apical end cap onto the cervix / PVW (tools/apex_close.py)
  L202_curvepress_apexclosed_P1  the same closure on L191_botfull_curvepress_cdamp1000_P1 (whole-curve pressure: the run
                                 that inverted from the apex): does the closed apex stop the inversion?
usage: py -3.10 build_batch202.py [K] [SRC NAME]"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from apex_close import close  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L192_botfull_distpress_cdamp1000_P1'
NAME = 'L202_distpress_apexclosed_P1'


def plane():
    ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot().find('Mesh')
    XR = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in ref.findall('Nodes') for n in b}
    P = np.array([XR[i] for i in range(23817, 24213)])
    c0 = P.mean(0)
    return c0, np.linalg.svd(P - c0)[2][-1]


if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    if len(sys.argv) > 3:
        SRC, NAME = sys.argv[2], sys.argv[3]
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    c0, nrm = plane()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    cx = [int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet346' for el in e for v in el.text.split(',')]
    nn = nrm if np.dot(np.mean([X[v] for v in cx], 0) - c0, nrm) < 0 else -nrm
    t, msg = close(t, K, c0, nn)
    print(msg)
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch202.py)\n{msg}\n')
