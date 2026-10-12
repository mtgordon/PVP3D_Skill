"""Batch 261 (2026-10-10, the user: "adopt your version as the tube main, and then run a version with pressure on distal half of
the curves as well as fully on the curves"). The tube main L258_tubecand_* has no pressure on the curves (removed for the
range of motion, L254: the curve pressure halves the apex motion). One change each, back on:
  full   Load-curves (0.014 MPa, lc 1, as L242_apexcontour_cx3_curvecx_attach_P1) on every Load-curves_surf face (436)
  dist   the same on the distal half only (face centroid height > -45 mm, 0 = hymen at rest; as build_batch192 B / L192)
NOT IN SOURCE (the curves are the tube's; their load as part of the canal's outer surface).
usage: py -3.10 build_batch261.py SRC NAME full|dist"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
CUT = -45.0
LOAD = ('\t\t<surface_load name="Load-curves" surface="Load-curves_surf" type="pressure">\n'
        '\t\t\t<pressure lc="1">0.014</pressure>\n\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n'
        '\t\t\t<linear>0</linear>\n\t\t\t<shell_bottom>0</shell_bottom>\n\t\t</surface_load>\n')

if __name__ == '__main__':
    src, name, mode = sys.argv[1:4]
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    assert '<surface_load name="Load-curves"' not in t
    a = t.index('\t\t<Surface name="Load-curves_surf">')
    b = t.index('</Surface>', a) + len('</Surface>\n')
    faces = [l for l in t[a:b].split('\n') if '<quad4' in l or '<tri3' in l]
    msg = f'NOT IN SOURCE: the pressure on the curves back on (Load-curves, 0.014 MPa) on all {len(faces)} Load-curves_surf faces'
    if mode == 'dist':
        mesh = ET.fromstring(t).find('Mesh')
        X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for bl in mesh.findall('Nodes') for n in bl}
        ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot().find('Mesh')
        XR = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for bl in ref.findall('Nodes') for n in bl}
        P = np.array([XR[i] for i in range(23817, 24213)])
        c0 = P.mean(0); nrm = np.linalg.svd(P - c0)[2][-1]
        cx = [int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet346' for el in e for v in el.text.split(',')]
        nn = nrm if np.dot(np.mean([X[v] for v in cx], 0) - c0, nrm) < 0 else -nrm
        keep = [l for l in faces
                if (np.mean([X[int(v)] for v in re.search(r'>([^<]+)<', l).group(1).split(',')], 0) - c0) @ nn > CUT]
        out = '\t\t<Surface name="Load-curves_surf">\n' + ''.join(
            re.sub(r'id="\d+"', f'id="{i}"', l) + '\n' for i, l in enumerate(keep, 1)) + '\t\t</Surface>\n'
        t = t[:a] + out + t[b:]
        msg = (f'NOT IN SOURCE: the pressure on the curves back on (Load-curves, 0.014 MPa) on the distal half only: {len(keep)} of '
               f'{len(faces)} faces (centroid height > {CUT:g} mm, 0 = hymen at rest)')
    i = t.index('\t\t<surface_load name="Load-PVW"')
    t = t[:i] + LOAD + t[i:]
    print(msg)
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch261.py)\n{msg}\n')
