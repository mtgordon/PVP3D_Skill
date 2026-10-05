"""Batch 118 (2026-09-29, the user's request): the canal as a tube, version 1: the tight fold (NOT IN SOURCE).
On the springs line's base L91_newline_vwyeoh_rhoi0_seamspr001 (the faithful walls + the 0.01 N/mm seam springs), the
AVW / cervix and the PVW are joined along both sides of the canal with shared nodes, the way the cervix shares nodes with
the PVW at the top:
  * every PVW lumen-edge node along the sides (tube.py: the PVW nodes nearest the 116 AVW / cervix seam nodes, 26-27 a side)
    is merged into its nearest AVW / cervix seam node (the PVW node moves onto it, <= ~2 mm; its references in elements,
    surfaces, node sets, spring sets and constraints are replaced); PVW nodes already shared with the cervix or the
    perineal body are left alone;
  * the AVW / cervix seam nodes between two merged ones (the AVW edge is twice as fine) follow the PVW's edge between them:
    linear constraints u - (1-s) u_A - s u_B = 0 per dof (s by arc length; penalty 10, maxaug 10 as LA_truss_ties);
  * the seam springs (constraint canal_seam_springs) removed: the fold replaces them.
Mechanically a line join along the lumen edge, like the seam tie of batch 90 (L86_newline_vwyeoh_rhoi0_seamtie), with the
gap closed and shared nodes.
  L118_tube_fold   L91_newline_vwyeoh_rhoi0_seamspr001, the seam springs -> the tight fold
usage: py -3.10 build_batch118.py"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch90 import seam_nodes, _bedges  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
LABEL = 'NOT IN SOURCE (the canal as a tube: tight fold)'


def replace_nodes(model, rep):
    """Replace node ids everywhere they are referenced: element and facet connectivity, node sets, spring sets and
    linear-constraint node attributes."""
    n = 0
    for el in model.root.iter():
        if el.tag in ('elem', 'delem', 'quad4', 'tri3', 'quad8', 'tri6') and el.text and ',' in el.text:
            ids = [int(v) for v in el.text.split(',')]
            if any(i in rep for i in ids):
                el.text = ','.join(str(rep.get(i, i)) for i in ids)
                n += 1
        elif el.tag == 'NodeSet' and el.text and el.text.strip():
            ids = [int(v) for v in re.split(r'[,\s]+', el.text.strip()) if v]
            if any(i in rep for i in ids):
                new, seen = [], set()
                for i in ids:
                    j = rep.get(i, i)
                    if j not in seen:
                        new.append(j)
                        seen.add(j)
                el.text = ','.join(map(str, new))
                n += 1
        elif el.tag in ('node', 'n') and el.get('id') and el.getparent() if hasattr(el, 'getparent') else False:
            pass
    for lc in model.root.iter('linear_constraint'):
        for nd in lc.findall('node'):
            i = int(nd.get('id'))
            if i in rep:
                nd.set('id', str(rep[i]))
                n += 1
    return n


def fold(model):
    path = model.src
    f = Feb(path)
    seam, fP, X = seam_nodes(Model7(path))
    eP = sorted({v for e in _bedges(fP) for v in e})
    PE = np.array([X[v] for v in eP])
    cx, peb = set(f.domain_nodes(tube.CX)), set(f.domain_nodes(tube.PEB))
    pvw_edge = sorted({eP[int(np.argmin(np.linalg.norm(PE - X[v], axis=1)))] for v in seam})
    pvw_edge = [p for p in pvw_edge if p not in cx and p not in peb]
    # merge each PVW edge node into its nearest seam node on the same side (one-to-one)
    rep, moved, used = {}, [], set()
    for p in sorted(pvw_edge, key=lambda p: min(np.linalg.norm(X[a] - X[p]) for a in seam)):
        side = [a for a in seam if np.sign(X[a][0]) == np.sign(X[p][0]) and a not in used]
        a = min(side, key=lambda a: np.linalg.norm(X[a] - X[p]))
        rep[p] = a
        used.add(a)
        moved.append(np.linalg.norm(X[a] - X[p]))
    # the seam springs out
    cons = model.root.find('Constraints')
    gone = [c.get('name') for c in cons.findall('constraint') if c.get('name') == 'canal_seam_springs']
    for c in list(cons):
        if c.get('name') == 'canal_seam_springs':
            cons.remove(c)
    nref = replace_nodes(model, rep)
    # the in-between seam nodes follow the merged ones along each side
    ties = []
    for sgn in (-1, 1):
        S = [a for a in seam if np.sign(X[a][0]) == sgn]
        # order along the side: the canal's long axis
        ax = np.linalg.svd(np.array([X[a] for a in S]) - np.mean([X[a] for a in S], axis=0))[2][0]
        S.sort(key=lambda a: X[a] @ ax)
        idx = [i for i, a in enumerate(S) if a in used]
        for i, a in enumerate(S):
            if a in used:
                continue
            lo = max([j for j in idx if j < i], default=None)
            hi = min([j for j in idx if j > i], default=None)
            if lo is None or hi is None:
                continue          # beyond the last merged node at a canal end: stays free (contact only)
            A, B = S[lo], S[hi]
            pts = S[lo:hi + 1]
            seg = np.r_[0, np.cumsum([np.linalg.norm(X[pts[k + 1]] - X[pts[k]]) for k in range(len(pts) - 1)])]
            s = seg[i - lo] / seg[-1]
            ties.append((a, A, 1 - s, B, s))
    c = ET.SubElement(cons, 'constraint', name='canal_fold_between', type='linear constraint')
    ET.SubElement(c, 'tol').text = '0.01'
    ET.SubElement(c, 'penalty').text = '10'
    ET.SubElement(c, 'maxaug').text = '10'
    for n, A, wa, B, wb in ties:
        for dof in ('x', 'y', 'z'):
            lc = ET.SubElement(c, 'linear_constraint')
            ET.SubElement(lc, 'node', id=str(n), bc=dof).text = '1'
            ET.SubElement(lc, 'node', id=str(A), bc=dof).text = f'{-wa:.6f}'
            ET.SubElement(lc, 'node', id=str(B), bc=dof).text = f'{-wb:.6f}'
    free = len(seam) - len(used) - len(ties)
    model.log.append(f'{LABEL}: the seam springs removed ({gone}); {len(rep)} PVW lumen-edge nodes along the sides merged '
                     f'into their nearest AVW / cervix seam nodes (moved median {np.median(moved):.2f}, max {max(moved):.2f} mm; '
                     f'{nref} references replaced; PVW nodes shared with the cervix / PeB left alone); {len(ties)} seam nodes '
                     f'between merged ones follow them (constraint canal_fold_between, penalty 10, maxaug 10); {free} seam '
                     f'nodes beyond the last merged node at a canal end left to the contact')
    return rep, ties


if __name__ == '__main__':
    name = 'L118_tube_fold'
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model7(os.path.join(RUNS, BASE, BASE + '.feb'))
    fold(m)
    emit(name, m)
