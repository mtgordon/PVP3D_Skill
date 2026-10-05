"""Batch 132 (2026-10-01, the user: "Let's keep 3 base versions when we make a new base version. 1 with just the springs, 1
with the springs and lofts (like 111), and one with lofts"). The current base L128_tube_r30s_pm (springs line: faithful
Yeoh walls, the 3 mm rounded tube, SlidingElastic1 seg_up 2, the PM as a structure with its anchors tied to it) is the
springs version. The other two:
  L132_chains_tube_pm   L128_tube_r30s_pm + every connector family as a chain of springs in a weak sheet (batch 111's
                        chains_all, sheet 5 % of E1, as L111_chains_all_m05). The PM anchors stay tied to the PM (each
                        chain starts at the old anchor node).
  L132_lofts_tube_pm    the lofts line with the faithful walls (L87_lofts_newline_vwyeoh_rhoi0: the user's four fitted
                        lofts AVW-Para, CL, USL, PM, the soft families as their connectors) + the same three changes as
                        the springs base: the 3 mm rounded tube (build_batch122.wrap, RIN 3; there are no seam springs on
                        this line), SlidingElastic1 seg_up 0 -> 2, and the PM as a structure (build_batch113: PM_Plane a
                        deformable 2 mm PM_Yeoh shell, outer arc clamped; the PM_PeB, PM_avw_bottom and PM-LA spring
                        anchors and the PM_fan loft's 26 origin nodes (BC-RP-PM-origins) tied to it, their BCs removed).
All NOT IN SOURCE as their parts are. usage: py -3.10 build_batch132.py [NAME ...]"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS, emit  # noqa: E402
import build_batch122 as b122  # noqa: E402
from build_batch113 import Model113, INNER, TIE, closest_on_tri, LABEL as PMLABEL  # noqa: E402
from build_batch113 import _pm_attach  # noqa: E402
from build_batch111 import Model111  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


class Model132(Model111, Model113):
    pass


def seg_up2(m):
    c = next(x for x in m.root.find('Contact') if x.get('name') == 'SlidingElastic1')
    assert c.findtext('seg_up').strip() == '0'
    c.find('seg_up').text = '2'
    m.log.append('solver-side: contact SlidingElastic1 seg_up 0 -> 2 (the tube\'s P1 speed fix, answer unchanged)')


def tie_nodeset_to_pm(m, ns):
    """tie every node of node set ns to the PM shell at its closest point (as build_batch113._pm_attach) and drop the
    set's fixed BC"""
    X = m.nodes()
    pm = [[int(v) for v in e.text.split(',')] for e in m.elem_blocks()['PM_Plane']]
    tris = []
    for c in pm:
        tris.append(c[:3])
        if len(c) == 4:
            tris.append([c[0], c[2], c[3]])
    A = np.array([X[t[0]] for t in tris]); B = np.array([X[t[1]] for t in tris]); C = np.array([X[t[2]] for t in tris])
    cen = (A + B + C) / 3
    ids = m.nodesets_ids(ns)
    c = next((x for x in m.root.find('Constraints') if x.get('name') == 'PM_anchor_ties'), None)
    assert c is not None
    off = []
    for n in ids:
        p = X[n]
        best = None
        for i in np.argsort(np.linalg.norm(cen - p, axis=1))[:30]:
            q, w = closest_on_tri(p, A[i], B[i], C[i])
            dd = np.linalg.norm(p - q)
            if best is None or dd < best[0]:
                best = (dd, i, w)
        dd, i, w = best
        off.append(dd)
        for dof in ('x', 'y', 'z'):
            lc = ET.SubElement(c, 'linear_constraint')
            ET.SubElement(lc, 'node', id=str(n), bc=dof).text = '1'
            for k in range(3):
                if w[k] > 1e-9:
                    ET.SubElement(lc, 'node', id=str(tris[i][k]), bc=dof).text = f'{-w[k]:.6f}'
    bnd = m.root.find('Boundary')
    for bc in list(bnd):
        if bc.get('node_set') == ns:
            bnd.remove(bc)
    m.log.append(f'{PMLABEL}: {ns} ({len(ids)} nodes, the PM_fan loft\'s anchor edge) tied to the PM shell (PM_anchor_ties); '
                 f'rest offset min / median / max {min(off):.2f} / {np.median(off):.2f} / {max(off):.2f} mm; its fixed BC '
                 f'removed')


def chains(m):
    m.chains_all(0.05)


def lofts(m):
    b122.RIN = 3.0
    _, jh, _ = b122.wrap(m)
    assert (jh > 0).all()
    seg_up2(m)
    m.pm_structure()
    _pm_attach(m, INNER + ('PM-LA-x%stiff_mat',))
    tie_nodeset_to_pm(m, 'BC-RP-PM-origins')


BUILDS = (('L132_chains_tube_pm', 'L128_tube_r30s_pm', chains),
          ('L132_lofts_tube_pm', 'L87_lofts_newline_vwyeoh_rhoi0', lofts))

if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        mdl = Model132(os.path.join(RUNS, base, base + '.feb'))
        fn(mdl)
        for line in mdl.log[-4:]:
            print('  ', line[:300])
        emit(name, mdl)
