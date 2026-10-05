"""Batch 141 (2026-10-01, the user: "it looks like the PVW part of the curved edges pushed through the AVW"). The lumen contact
SlidingElastic1 covers only the walls' lumen faces (PVW: primary; AVW / cervix: secondary); the curved edges' inner
(channel) faces are in no contact, so when the channel closes and the curves fold, the PVW half of a curve passes through
the AVW and the AVW half of the curve. Fix (NOT IN SOURCE, as the curves are): the curves' inner faces join the lumen
contact: the PVW half's (canal_wrap_hex_pvw) inner faces are added to SlidingElastic1Primary, the AVW half's
(canal_wrap_hex_avw) to SlidingElastic1Secondary. Inner faces = the wrap's boundary faces whose outward normal points
toward the local centre of the U (the midpoint of the AVW seam node and its paired PVW lumen-edge node) within 60 degrees;
outward winding, as the walls' lumen faces.
  L141_thin275_narrow3_wrapcontact         L137_tube_thin275_narrow3_flare_pmin + the curves in the lumen contact
  L141_thin275_narrow3_wrapP_wrapcontact   L138_tube_thin275_narrow3_flare_pmin_wrapP (+ the pressure on the curves) + the same
usage: py -3.10 build_batch141.py"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS  # noqa: E402
from build_batch90 import seam_nodes  # noqa: E402
from build_batch118 import fold  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
ORIG = os.path.join(RUNS, 'L91_newline_vwyeoh_rhoi0_seamspr001', 'L91_newline_vwyeoh_rhoi0_seamspr001.feb')
BUILDS = (('L141_thin275_narrow3_wrapcontact', 'L137_tube_thin275_narrow3_flare_pmin'),
          ('L141_thin275_narrow3_wrapP_wrapcontact', 'L138_tube_thin275_narrow3_flare_pmin_wrapP'))


def inner_faces(f, dom, C):
    X = f.nodes
    walln = {n for dm in (tube.AVW, tube.CX, tube.PVW) for n in f.domain_nodes(dm)}
    out = []
    for fc in tube.boundary_faces(f, (dom,)):
        if all(n in walln for n in fc):
            continue
        P = np.array([X[n] for n in fc])
        cen = P.mean(0)
        nrm = np.cross(P[2] - P[0], P[3] - P[1])
        nrm /= np.linalg.norm(nrm)
        c = C[int(np.argmin(np.linalg.norm(C - cen, axis=1)))]
        r = (c - cen) / np.linalg.norm(c - cen)
        if nrm @ r > np.cos(np.radians(60)):
            out.append(fc)
    return out


def add_faces(t, surf, faces):
    i = t.index(f'<Surface name="{surf}">')
    j = t.index('</Surface>', i)
    blk = t[i:j]
    ids = [int(v) for v in re.findall(r'<(?:quad4|tri3) id="(\d+)"', blk)]
    k0 = max(ids) + 1
    add = ''.join(f'\t\t\t<quad4 id="{k0 + q}">{",".join(map(str, fc))}</quad4>\n' for q, fc in enumerate(faces))
    return t[:j] + add.lstrip('\t') + '\t\t' + t[j:] if blk.endswith('\n') else t[:j] + '\n' + add + '\t\t' + t[j:]


if __name__ == '__main__':
    seam, _, _ = seam_nodes(Model7(ORIG))
    rep, _ = fold(Model7(ORIG))
    inv = {a: p for p, a in rep.items()}
    for name, base in BUILDS:
        d = os.path.join(RUNS, name)
        assert not os.path.exists(d), f'{name} exists; not overwriting'
        src = os.path.join(RUNS, base, base + '.feb')
        f, t = Feb(src), open(src, encoding='utf-8').read()
        C = np.array([0.5 * (np.array(f.nodes[a]) + np.array(f.nodes[inv[a]])) for a in seam if a in inv])
        fp = inner_faces(f, 'canal_wrap_hex_pvw', C)
        fa = inner_faces(f, 'canal_wrap_hex_avw', C)
        t2 = add_faces(add_faces(t, 'SlidingElastic1Primary', fp), 'SlidingElastic1Secondary', fa)
        os.makedirs(d)
        open(os.path.join(d, name + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
        msg = (f'{name}: {base} + NOT IN SOURCE: the curved edges\' inner (channel) faces in the lumen contact '
               f'SlidingElastic1: {len(fp)} PVW-half faces added to SlidingElastic1Primary, {len(fa)} AVW-half faces to '
               f'SlidingElastic1Secondary\n')
        open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
        print(msg)
