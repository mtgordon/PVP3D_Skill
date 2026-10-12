"""Batch 192 (2026-10-06, the user: "yes, start both"). L191 (equal PM-arch springs + curve pressure + curve damping C 1000)
inverted from the APICAL end (the user, in FEBio Studio) at t ~0.78 and blew up, instead of pushing out / turning inside
out at the distal end. Two single-change diagnostics from L191, P1:
  L192_curvepress_cdamp1000_P1        A: the PM-arch springs back to L169's (PM_avw_bottom at the source law x1e-5, i.e.
                                         L189 + curve damping C 1000): are the stiff distal anchors what stop the eversion?
                                         (NOT IN SOURCE; breaks the equal-arch-spring rule ON PURPOSE, as a diagnostic)
  L192_botfull_distpress_cdamp1000_P1 B: the curve pressure only on faces with centroid height > -45 mm (the distal half;
                                         0 = hymen at rest): does loading the apical curve faces drive the inversion?
usage: py -3.10 build_batch192.py"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L191_botfull_curvepress_cdamp1000_P1'
CUT = -45.0


def damp1000(t):
    for d in ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw'):
        a = t.index(f'<body_load name="wall_damping_{d}"')
        b = t.index('</body_load>', a)
        blk, k = re.subn(r'<C>200</C>', '<C>1000</C>', t[a:b])
        assert k == 1, d
        t = t[:a] + blk + t[b:]
    return t


if __name__ == '__main__':
    # A: L189 (L169 arch springs + curve pressure) + curve damping 1000
    t = damp1000(open(os.path.join(RUNS, 'L189_selfall_curvepress_P1', 'L189_selfall_curvepress_P1.feb'), encoding='utf-8').read())
    write('L192_curvepress_cdamp1000_P1', t,
          'base: claude_diag/runs/L189_selfall_curvepress_P1/L189_selfall_curvepress_P1.feb (tools/build_batch192.py)\n'
          'NOT IN SOURCE: the curves\' mass damping C 200 -> 1000 /s; = L191 with the PM-arch springs as L169 '
          '(PM_avw_bottom at the source law x1e-5): a diagnostic, breaks the equal-arch-spring rule on purpose\n')
    # B: L191 with the curve pressure on the distal half only
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot().find('Mesh')
    XR = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in ref.findall('Nodes') for n in b}
    P = np.array([XR[i] for i in range(23817, 24213)])
    c0 = P.mean(0)
    nrm = np.linalg.svd(P - c0)[2][-1]
    cx = [int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet346' for el in e for v in el.text.split(',')]
    nn = nrm if np.dot(np.mean([X[v] for v in cx], 0) - c0, nrm) < 0 else -nrm
    a = t.index('\t\t<Surface name="Load-curves_surf">')
    b = t.index('</Surface>', a) + len('</Surface>\n')
    lines = t[a:b].split('\n')
    faces = [l for l in lines if '<quad4' in l or '<tri3' in l]
    keep = [l for l in faces
            if (np.mean([X[int(v)] for v in re.search(r'>([^<]+)<', l).group(1).split(',')], 0) - c0) @ nn > CUT]
    out = '\t\t<Surface name="Load-curves_surf">\n' + ''.join(
        re.sub(r'id="\d+"', f'id="{i}"', l) + '\n' for i, l in enumerate(keep, 1)) + '\t\t</Surface>\n'
    t = t[:a] + out + t[b:]
    write('L192_botfull_distpress_cdamp1000_P1', t,
          f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch192.py)\n'
          f'NOT IN SOURCE: Load-curves on the distal half only: {len(keep)} of {len(faces)} faces (centroid height > {CUT} mm, '
          f'0 = hymen at rest)\n')
    print('B keeps', len(keep), 'of', len(faces))
