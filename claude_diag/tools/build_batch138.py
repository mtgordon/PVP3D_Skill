"""Batch 138 (2026-10-01, the user: "add the pressure to the outside of the curved edges surface as well and run it").
L137_tube_thin275_narrow3_flare_pmin + NOT IN SOURCE pressure Load-wrap: the source's wall pressure (0.014 MPa, load curve 1,
as Load-AVW / Load-PVW, which sit on the walls' outer surfaces with their normals pointing away from the lumen) on the
curved edges' outer (convex) faces: the boundary faces of canal_wrap_hex_avw / _pvw whose outward normal points away
from the local centre of the U (the midpoint of the AVW seam node and its paired PVW lumen-edge node, nearest the face)
within 60 degrees of the radial direction; the inner (channel) faces and the end faces at the first / last station are
left out. Facets written with outward winding, as the walls' load surfaces.
  L138_tube_thin275_narrow3_flare_pmin_wrapP
usage: py -3.10 build_batch138.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS  # noqa: E402
from build_batch90 import seam_nodes  # noqa: E402
from build_batch118 import fold  # noqa: E402
from build_batch131 import facets_xml  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE, NAME = 'L137_tube_thin275_narrow3_flare_pmin', 'L138_tube_thin275_narrow3_flare_pmin_wrapP'
ORIG = os.path.join(RUNS, 'L91_newline_vwyeoh_rhoi0_seamspr001', 'L91_newline_vwyeoh_rhoi0_seamspr001.feb')

if __name__ == '__main__':
    src = os.path.join(RUNS, BASE, BASE + '.feb')
    d = os.path.join(RUNS, NAME)
    assert not os.path.exists(d), f'{NAME} exists; not overwriting'
    f, t = Feb(src), open(src, encoding='utf-8').read()
    X = {n: np.array(p, float) for n, p in f.nodes.items()}
    seam, _, _ = seam_nodes(Model7(ORIG))
    rep, _ = fold(Model7(ORIG))
    inv = {a: p for p, a in rep.items()}
    C = np.array([0.5 * (X[a] + X[inv[a]]) for a in seam if a in inv])        # local centres of the U
    faces = tube.boundary_faces(f, ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw'))
    walln = {n for dm in (tube.AVW, tube.CX, tube.PVW) for n in f.domain_nodes(dm)}
    out = []
    for fc in faces:
        if all(n in walln for n in fc):          # the U's ends lying on the walls' side faces
            continue
        P = np.array([X[n] for n in fc])
        cen = P.mean(0)
        nrm = np.cross(P[2] - P[0], P[3] - P[1])
        nrm /= np.linalg.norm(nrm)
        c = C[int(np.argmin(np.linalg.norm(C - cen, axis=1)))]
        r = cen - c
        r /= np.linalg.norm(r)
        if nrm @ r > np.cos(np.radians(60)):
            out.append(fc)
    area = sum(0.5 * np.linalg.norm(np.cross(X[fc[2]] - X[fc[0]], X[fc[3]] - X[fc[1]])) for fc in out)
    i = t.index('<Surface name="Load-AVW_surf">')
    t2 = t[:i] + facets_xml('Load-wrap_surf', out).lstrip('\t') + '\t\t' + t[i:]
    j = t2.index('<surface_load name="Load-AVW"')
    k = t2.index('</surface_load>', j) + len('</surface_load>')
    load = t2[j:k].replace('name="Load-AVW" surface="Load-AVW_surf"', 'name="Load-wrap" surface="Load-wrap_surf"')
    assert 'Load-wrap_surf' in load
    t2 = t2[:k] + '\n\t\t' + load + t2[k:]
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
    msg = (f'{NAME}: {BASE} + NOT IN SOURCE: pressure Load-wrap (0.014 MPa, lc 1, as Load-AVW) on the curved edges\' outer '
           f'faces: {len(out)} of the {len(faces)} wrap boundary faces, {area:.0f} mm2\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
