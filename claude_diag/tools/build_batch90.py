"""Batch 90 (2026-09-27 night): the canal seam tie (NOT IN SOURCE; for the user to adopt or not), on the batch-88 leader.

The canal's contact (SlidingElastic1: PVW front + PeB primary, AVW + cervix secondary) closes the canal at its sides by
contact alone: the AVW and PVW share no nodes, and the secondary's lateral edges lie 0.12-0.36 mm from the primary's edge
along the whole canal (the seam; tools/build_batch79.py seam_out finds 145 edge nodes within 0.75 mm: both lateral edges,
cervix to introitus, and the cervix's proximal end row). The faithful (Yeoh) walls crawl where these edge nodes flip on and
off the other surface's last facet. The fix ties each lateral-seam secondary node to the nearest point of the primary
surface at its rest offset, with linear constraints u_node - sum_k w_k u_k = 0 per dof (w_k: the facet's bilinear shape
functions at that point), as tested on claude_diag/seam_tie/seam_mini.py (offset kept to 1e-5 mm; penalty 10 as the model's
LA_truss_ties). The proximal end row (the cervix, already joined to the PVW by shared nodes) is left as it is.
  L86_newline_vwyeoh_rhoi0_seamtie   L85_newline_vwyeoh_rhoi0 + the lateral seam tie
usage: py -3.10 build_batch90.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from variants7 import Model7, JOBS, emit
from build_batch79 import _bedges, _dseg

sys.path.insert(0, os.path.join(JOBS, 'claude_diag', 'seam_tie'))
from seam_mini import project  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def seam_nodes(m, tol=0.75, pair='SlidingElastic1', lateral=True):
    """secondary edge nodes within tol of the primary's edge; lateral=True drops the canal's two end rows (corners kept)."""
    X = m.nodes()
    surf = {s.get('name'): s for s in m.mesh.findall('Surface')}
    fl = lambda s: [tuple(int(v) for v in f.text.split(',')) for f in s]
    fS, fP = fl(surf[pair + 'Secondary']), fl(surf[pair + 'Primary'])
    eP = _bedges(fP)
    bS = sorted({v for e in _bedges(fS) for v in e})
    seam = [v for v in bS if min(_dseg(X[v], X[a], X[b]) for a, b in eP) < tol]
    if lateral:
        C = np.array([X[v] for v in sorted({v for f in fS for v in f})])
        c0 = C.mean(0)
        ax = np.linalg.svd(C - c0, full_matrices=False)[2]
        p = {v: (X[v] - c0) @ ax.T for v in bS}
        lo, hi = min(q[0] for q in p.values()), max(q[0] for q in p.values())
        def keep(v):
            if lo + 0.5 < p[v][0] < hi - 0.5:
                return True
            row = [w for w in bS if abs(p[w][0] - p[v][0]) < 0.5]          # an end row: keep only its two corners
            return abs(p[v][1]) >= max(abs(p[w][1]) for w in row) - 1e-6
        seam = [v for v in seam if keep(v)]
    return seam, fP, X


def seam_tie(m, penalty=10, maxaug=10, **kw):
    """maxaug 10: a tie (augmented Lagrangian, offset kept to ~1e-5 mm). maxaug 0: no augmentation, so each constraint is
    a zero-length spring of stiffness `penalty` N/mm from the node to the facet point, in every direction (the seam
    springs; checked on claude_diag/seam_tie/seam_mini.py: the stretch grows as the stiffness drops)."""
    seam, fP, X = seam_nodes(m, **kw)
    cen = np.array([np.mean([X[n] for n in f], axis=0) for f in fP])
    ties, offs = [], []
    for v in seam:
        near = [fP[i] for i in np.argsort(np.linalg.norm(cen - X[v], axis=1))[:12]]
        d, f, w = project(X[v], X, near)
        ties.append((v, f, w)); offs.append(d)
    if maxaug == 0:
        return _seam_constraints(m, ties, offs, penalty, 0, 'canal_seam_springs',
                                 f'NOT IN SOURCE (canal seam springs): {len(ties)} SlidingElastic1Secondary lateral-edge '
                                 f'nodes (AVW / cervix) each joined to the point of SlidingElastic1Primary (PVW) it faces at '
                                 f'rest ({min(offs):.2f}-{max(offs):.2f} mm, median {np.median(offs):.2f}) by a zero-length '
                                 f'spring of {penalty} N/mm in every direction (penalty-only linear constraints, maxaug 0)')
    return _seam_constraints(m, ties, offs, penalty, maxaug, 'canal_seam_tie',
                             f'NOT IN SOURCE (canal seam tie): {len(ties)} SlidingElastic1Secondary lateral-edge nodes (AVW / '
                             f'cervix) tied to the nearest point of SlidingElastic1Primary (PVW / PeB) at their rest offset '
                             f'({min(offs):.2f}-{max(offs):.2f} mm, median {np.median(offs):.2f}) by linear constraints '
                             f'(penalty {penalty}, tol 0.01, maxaug {maxaug})')


def _seam_constraints(m, ties, offs, penalty, maxaug, name, note):
    cons = m.root.find('Constraints')
    c = ET.SubElement(cons, 'constraint', {'name': name, 'type': 'linear constraint'})
    for k, val in (('tol', 0.01), ('penalty', penalty), ('maxaug', maxaug)):
        ET.SubElement(c, k).text = str(val)
    for v, f, w in ties:
        for dof in 'xyz':
            lc = ET.SubElement(c, 'linear_constraint')
            ET.SubElement(lc, 'node', {'id': str(v), 'bc': dof}).text = '1'
            for n, wk in zip(f, w):
                if abs(wk) > 1e-9:
                    ET.SubElement(lc, 'node', {'id': str(n), 'bc': dof}).text = '%.9g' % -wk
    m.log.append(note)
    return len(ties)


BUILDS = (('L86_newline_vwyeoh_rhoi0_seamtie', 'L85_newline_vwyeoh_rhoi0', seam_tie),)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        print(m.log[-1])
        emit(name, m)
