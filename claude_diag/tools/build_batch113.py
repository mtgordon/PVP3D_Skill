"""Batch 113 (2026-09-29, the user's request): the perineal membrane as a structure instead of a reference.
Our model already has the PM: PM_Plane (the source's display body, a U-shaped band of 350 S4R / quad4 shells; in FEBio a
fixed rigid body, 0.5 mm, used only as the hymenal reference). The user pointed to OPAL325_PM_mid in
Normal_Generic_copy.inp (a different model; NOT brought in) for how a PM is handled there:
  * its outer edge is fixed: *Boundary PM on nodes 1-53, dofs 1-6 (clamped), one continuous 53-node stretch of its
    121-node boundary; the perineal body is tied to 5 consecutive nodes (101-105) of the free part (Tie PM-PBody);
  * shell section 2.0 mm, material PARAVAG_H_highdensity: *Hyperelastic, n=2 (polynomial) from uniaxial test data,
    0.0027 MPa at 20 % ... 0.031 MPa at 100 % nominal strain (the 10 % point, 0.0269, is 10x the 20 % one: a likely
    typo; the fit is the same with or without it), density 9e-07.
Step 1 (this batch): PM_Plane becomes a deformable shell, clamped along its outer arc, with nothing attached:
  * material PM_Yeoh: uncoupled Yeoh c1 0.002770494, c2 0.001810007 (skill scripts/fit_yeoh.py to the test data without
    the 10 % point, N=2, max error 19 %; an I1-only stand-in for Abaqus' polynomial N=2, a fidelity approximation),
    k 0.25 (Abaqus' constant-nu volume change at ~50 % strain; k / mu0 = 45), density 1.06e-9 (tissue; the source's
    9e-07 is a mass-scaled "highdensity" value), thickness 2.0 mm, shell_normal_nodal 0;
  * the outer arc: the boundary loop's stretch between its two outer corners (the four sharpest corners split the U into
    outer arc, inner arc and two short ends), the one farther from the vaginal walls and perineal body; x, y, z and the
    back face's sx, sy, sz fixed (clamped, as dofs 1-6);
  * the rigid-body constraint on PM_Plane removed. The POP-Q tools keep taking the hymenal plane from PM_Plane's rest
    coordinates (pop_q.py), so the measuring plane does not move.
Nothing attaches to it yet, so the answer should equal the springs line's: the test of the conversion.
  L113_springs_pm_step1   L87_springs_newline_rhoi0 + the PM as a clamped deformable shell (nothing attached)
usage: py -3.10 build_batch113.py [NAME ...]
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from febmodel import Feb  # noqa: E402
import pm_survey  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L87_springs_newline_rhoi0'
WANT = set(sys.argv[1:])
LABEL = 'NOT IN SOURCE (the PM as a structure)'
PM_YEOH = dict(c1='0.002770494', c2='0.001810007', k='0.25', density='1.06e-09')
T_PM = '2'


def pm_arcs(path):
    """(loop, outer arc node ids, inner arc node ids) of PM_Plane in the model at path."""
    f = Feb(path)
    loop, tris = pm_survey.pm_geometry(f)
    P = np.array([f.nodes[n] for n in loop])
    L = len(loop)
    v1 = np.roll(P, 1, 0) - P
    v2 = np.roll(P, -1, 0) - P
    ang = np.degrees(np.arccos(np.clip((v1 * v2).sum(1) / np.linalg.norm(v1, axis=1) / np.linalg.norm(v2, axis=1), -1, 1)))
    corners = sorted(np.argsort(ang)[:4])
    tis = np.array([f.nodes[n] for dom in ('_PickedSet347', '_PickedSet64', '_PickedSet66') for n in f.domain_nodes(dom)])
    dt = np.array([np.linalg.norm(tis - p, axis=1).min() for p in P])
    stretches = [[(c0 + s) % L for s in range((c1 - c0) % L + 1)] for c0, c1 in zip(corners, corners[1:] + corners[:1])]
    long2 = sorted(stretches, key=len)[-2:]              # the two arcs (the ends are short)
    outer = max(long2, key=lambda s: dt[s].mean())
    inner = min(long2, key=lambda s: dt[s].mean())
    return loop, [loop[i] for i in outer], [loop[i] for i in inner], corners


class Model113(Model7):
    def pm_structure(self):
        path = self.src
        loop, outer, inner, corners = pm_arcs(path)
        mats = self.root.find('Material')
        old = next(m for m in mats if m.get('name') == 'PM_Plane_display')
        mid = old.get('id')
        i = list(mats).index(old)
        mats.remove(old)
        m = ET.Element('material', id=mid, name='PM_Yeoh', type='Yeoh')
        for k in ('density', 'c1', 'c2', 'k'):
            ET.SubElement(m, k).text = PM_YEOH[k]
        mats.insert(i, m)
        dom = next(d for d in self.root.find('MeshDomains') if d.get('name') == 'PM_Plane')
        dom.set('mat', 'PM_Yeoh')
        dom.find('shell_thickness').text = T_PM
        if dom.find('shell_normal_nodal') is None:
            ET.SubElement(dom, 'shell_normal_nodal').text = '0'
        # the rigid-body constraint on the display body
        gone = []
        rig = self.root.find('Rigid')
        if rig is not None:
            for r in list(rig):
                if r.findtext('rb') in ('PM_Plane_display', mid):
                    rig.remove(r)
                    gone.append(r.get('name') or r.tag)
            if len(rig) == 0:
                self.root.remove(rig)
        # the outer arc clamped
        ns = ET.Element('NodeSet', name='BC-PM-outer-arc')
        ns.text = ','.join(map(str, outer))
        last_ns = self.mesh.findall('NodeSet')[-1]
        self.mesh.insert(list(self.mesh).index(last_ns) + 1, ns)
        bnd = self.root.find('Boundary')
        bc = ET.SubElement(bnd, 'bc', name='BC-PM-outer-arc', node_set='BC-PM-outer-arc', type='zero displacement')
        for d in ('x_dof', 'y_dof', 'z_dof'):
            ET.SubElement(bc, d).text = '1'
        bc2 = ET.SubElement(bnd, 'bc', name='BC-PM-outer-arc-shell', node_set='BC-PM-outer-arc', type='zero shell displacement')
        for d in ('sx_dof', 'sy_dof', 'sz_dof'):
            ET.SubElement(bc2, d).text = '1'
        self.log.append(f'{LABEL}: PM_Plane from the fixed rigid display body (material PM_Plane_display) to a deformable '
                        f'shell: material PM_Yeoh (Yeoh c1 {PM_YEOH["c1"]}, c2 {PM_YEOH["c2"]}, k {PM_YEOH["k"]}, density '
                        f'{PM_YEOH["density"]}; fitted to OPAL325_PM_mid\'s PARAVAG_H_highdensity uniaxial data), thickness '
                        f'{T_PM} mm (PM_mid\'s section), shell_normal_nodal 0; rigid constraints removed: {gone or "none"}; '
                        f'outer arc clamped (x, y, z, sx, sy, sz): {len(outer)} of the {len(loop)} boundary nodes (corners at '
                        f'loop positions {corners}), node set BC-PM-outer-arc; inner arc {len(inner)} nodes free')
        return outer, inner


def closest_on_tri(p, a, b, c):
    """Closest point to p on triangle abc and its barycentric weights (wa, wb, wc)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = ab @ ap, ac @ ap
    if d1 <= 0 and d2 <= 0:
        return a, (1, 0, 0)
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return b, (0, 1, 0)
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return a + v * ab, (1 - v, v, 0)
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return c, (0, 0, 1)
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return a + w * ac, (1 - w, 0, w)
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b), (0, 1 - w, w)
    den = 1 / (va + vb + vc)
    v, w = vb * den, vc * den
    return a + ab * v + ac * w, (1 - v - w, v, w)


# spring sets whose anchor (the end in a fixed BC-RP-PM* set) is tied to the PM in step 2 (the user's choice, ~08:50)
INNER = ('PM_PeB_Left_conn', 'PM_PeB_Right_conn', 'PM_avw_bottom_left_conn', 'PM_avw_bottom_right_conn')
OUTER_BOTTOM = ('PM_conn', 'PM-LA-x%stiff_mat')
TIE = {'penalty': '100', 'maxaug': '0'}


def _pm_attach(model, sets):
    """Each spring set's anchor node (the end in a fully fixed BC set) loses its BC and follows the PM shell at its
    closest point: per dof, u_anchor - sum_k w_k u_k = 0 over the PM facet's nodes (weights = barycentric, a quad split
    in two triangles), so the anchor keeps its rest offset (0-2 mm); penalty 100 N/mm, maxaug 0 (a stiff zero-length
    spring, skill gotcha 27). Orphan nodes of the emptied BC sets (in no element or spring) stay fixed."""
    X = model.nodes()
    pm = [[int(v) for v in e.text.split(',')] for e in model.elem_blocks()['PM_Plane']]
    tris = []
    for c in pm:
        tris.append(c[:3])
        if len(c) == 4:
            tris.append([c[0], c[2], c[3]])
    A = np.array([X[t[0]] for t in tris]); B = np.array([X[t[1]] for t in tris]); C = np.array([X[t[2]] for t in tris])
    cen = (A + B + C) / 3
    bcset = {}
    for bc in model.root.find('Boundary'):
        ns = bc.get('node_set')
        if ns and ns.startswith('BC-RP-PM') and bc.get('type') == 'zero displacement':
            for n in model.nodesets_ids(ns):
                bcset[n] = ns
    used = set()
    for blk in model.elem_blocks().values():
        for e in blk:
            used.update(int(v) for v in e.text.split(','))
    for ds in model.mesh.findall('DiscreteSet'):
        for e in ds.findall('delem'):
            used.update(int(v) for v in e.text.split(','))
    anchors = {}
    for s in sets:
        ds = next(d for d in model.mesh.findall('DiscreteSet') if d.get('name') == s)
        for e in ds.findall('delem'):
            a, b = (int(v) for v in e.text.split(','))
            anc = [n for n in (a, b) if n in bcset]
            assert len(anc) == 1, (s, a, b)
            anchors[anc[0]] = s
    cons, off = [], []
    for n in sorted(anchors):
        p = X[n]
        best = None
        for i in np.argsort(np.linalg.norm(cen - p, axis=1))[:30]:
            q, w = closest_on_tri(p, A[i], B[i], C[i])
            dd = np.linalg.norm(p - q)
            if best is None or dd < best[0]:
                best = (dd, i, w)
        dd, i, w = best
        cons.append((n, [(tris[i][k], w[k]) for k in range(3) if w[k] > 1e-9]))
        off.append(dd)
    c = ET.SubElement(model.root.find('Constraints'), 'constraint', name='PM_anchor_ties', type='linear constraint')
    ET.SubElement(c, 'tol').text = '0.01'
    ET.SubElement(c, 'penalty').text = TIE['penalty']
    ET.SubElement(c, 'maxaug').text = TIE['maxaug']
    for n, ws in cons:
        for dof in ('x', 'y', 'z'):
            lc = ET.SubElement(c, 'linear_constraint')
            ET.SubElement(lc, 'node', id=str(n), bc=dof).text = '1'
            for m_, wk in ws:
                ET.SubElement(lc, 'node', id=str(m_), bc=dof).text = f'{-wk:.6f}'
    # the anchors out of their BC sets; a set left with only orphan nodes keeps them fixed
    changed = []
    for ns in sorted(set(bcset[n] for n in anchors)):
        ids = model.nodesets_ids(ns)
        keep = [n for n in ids if n not in anchors]
        node = model._nodeset(ns)
        for ch in list(node):
            node.remove(ch)
        node.text = ','.join(map(str, keep))
        orphans = [n for n in keep if n not in used]
        if not keep:
            bnd = model.root.find('Boundary')
            for bc in list(bnd):
                if bc.get('node_set') == ns:
                    bnd.remove(bc)
            model.mesh.remove(node)
        changed.append(f'{ns}: {len(ids)} -> {len(keep)} nodes ({len(orphans)} of them in no element or spring)')
    per = {}
    for n, s in anchors.items():
        per.setdefault(s, []).append(n)
    model.log.append(f'{LABEL}: anchors tied to the PM shell (constraint PM_anchor_ties, penalty {TIE["penalty"]} N/mm, '
                     f'maxaug {TIE["maxaug"]}): ' + ', '.join(f'{s} {len(v)}' for s, v in per.items())
                     + f'; rest offset from the PM surface min / median / max {min(off):.2f} / {np.median(off):.2f} / '
                     f'{max(off):.2f} mm; fixed BCs removed from them: {"; ".join(changed)}')


def step2(sets):
    def fn(m):
        m.pm_structure()
        _pm_attach(m, sets)
    return fn


BUILDS = [('L113_springs_pm_step1', BASE, lambda m: m.pm_structure()),
          # step 2 (the user's choice ~08:50: the inner-arc and the outer-bottom anchors; tissue density)
          ('L113_springs_pm_step2', BASE, step2(INNER + OUTER_BOTTOM)),
          ('L113_springs_pm_step2_inner', BASE, step2(INNER)),
          # step 2 on the paper's cases (the user, ~09:50: "get the testing started on the PM changes")
          ('L113_springs_pm_step2_paperP1', 'L107_springs_newline_rhoi0_paperP1', step2(INNER + OUTER_BOTTOM)),
          ('L113_springs_pm_step2_paperP2', 'L107_springs_newline_rhoi0_paperP2', step2(INNER + OUTER_BOTTOM))]

if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        if os.path.exists(os.path.join(RUNS, name)):
            print(f'{name} exists; not overwriting')
            continue
        mdl = Model113(os.path.join(RUNS, base, base + '.feb'))
        fn(mdl)
        emit(name, mdl)
