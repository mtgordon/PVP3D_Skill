"""Batch 109 (2026-09-29, the user's request of 2026-09-28): the third line, "fibre lofts".
Every family ever lofted (the 14 *_fan shell domains of L9_fitall_edge) back on the springs line L87_springs_newline_rhoi0
in place of its springs, as a weak matrix + fibres along the source connector lines that carry the family's connector
law (tools/fibre_lofts.py):
  * elements, ShellDomain (thickness 0.49, shell_normal_nodal 0) and the loft ties (AVW-Para-L, USL-L tied-node-on-facet)
    from L9_fitall_edge (the same nodes and ids as the springs line; the loft nodes shared with the tissue, the arcus
    chain and the BC sets are the same in both); the family's springs removed; the CL / USL BC sets back to the lofts
    model's (Model7.connectors_to_lofts);
  * material "<loft>_fib": solid mixture (density 1.06e-9, tissue) of a neo-Hookean matrix (E = m x E1, v 0.3; floor
    m x 0.107 MPa, the P-arcus E1, for the near-zero families and, in the paper's cases, the impaired ones) + fiber-pow-linear (E1, beta 2, lam0 1.005: P ~ E1 (lam-1),
    the tables' first segment) + fiber-exp-pow (ksi2, alpha, beta, lam0 1: the stiffening), fitted by the strip model to
    the family's connector tables over u up to ~1.1x the paper case P2's largest elongation;
  * fibre directions per element: a mat_axis ElementData (a = the blend of the two nearest connector lines' directions,
    projected into the element's plane; d = the element normal); the fibres' <fiber type="vector">1,0,0</fiber> reads it.
Names: L109_fibre_<family>_m<matrix %>; "all" = the 14 lofts. NOT IN SOURCE: the lofts themselves (the source has
connectors); the laws are fitted to the source's connector tables.
usage: py -3.10 build_batch109.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import copy
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit, _insert_like  # noqa: E402
from variants6 import Model6, LINE_BC  # noqa: E402
from variants2 import renumber_discrete_materials  # noqa: E402
import fibre_lofts as fl  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L87_springs_newline_rhoi0'
WANT = set(sys.argv[1:])
LABEL = 'NOT IN SOURCE (fibre loft in place of the connectors)'
# family group -> (lofts, connector families as lofts_to_connectors named them, spring set name patterns)
GROUPS = {
    'avwpara': (('AVW-Para-L_fan', 'AVW-Para-R_fan'), ('AVW-Para-L', 'AVW-Para-R')),
    'cl': (('CL-L_fan', 'CL-R_fan'), ('CL-L', 'CL-R')),
    'usl': (('USL-L_fan', 'USL-R_fan'), ('USL-L', 'USL-R')),
    'pm': (('PM_fan',), ('PM',)),
    'parcus': (('P-arcus-L_fan', 'P-arcus-R_fan'), ('Parcus',)),
    'pmpeb': (('PM_PeB_Left_fan', 'PM_PeB_Right_fan'), ('PM_PeB_Left_', 'PM_PeB_Right_')),
    'avwbot': (('PM_avw_bottom_left_fan', 'PM_avw_bottom_right_fan'), ('PM_avw_bottom_left_', 'PM_avw_bottom_right_')),
    'pebconstrin': (('PeB-constrin_fan',), ('PeB-constrin',)),
}
# fit range: ~1.1x the largest connector elongation of the paper's case P2 on the springs line
# (L107_springs_newline_rhoi0_paperP2), at least 10 mm
UHI = {'AVW-Para': 10, 'CL': 42, 'USL': 39, 'PM_fan': 10, 'PM_PeB': 14, 'PM_avw': 32, 'PeB-c': 10, 'P-arcus': 35}
E1_FLOOR = 0.107          # the P-arcus E1: the matrix floor is m x this for the near-zero families

_FAMS = None


def fams():
    global _FAMS
    if _FAMS is None:
        _FAMS = fl.families()
    return _FAMS


def law(loft, scale=1.0):
    """(E1, ksi2, alpha, beta, rms) of the loft's two-fibre fit; E1 and ksi2 x scale (a paper case's impairment)."""
    key = next(k for k in UHI if loft.startswith(k))
    u = np.linspace(UHI[key] / 100, UHI[key], 100)
    T, (err, E1, k2, a, b) = fl.fit2(fams()[loft], u)
    return E1 * scale, k2 * scale, a, b, err


def lam_cap(loft):
    """The largest stretch in the loft's fit range: 1 + u_max / (its shortest connector)."""
    key = next(k for k in UHI if loft.startswith(k))
    return 1 + UHI[key] / float(np.min(fams()[loft]['L']))


CAP = {'on': False}      # batch 117: the stiffening fibre linear beyond the fit range (fiber-exp-pow-linear)


def material_xml(name, mid, E1, ksi2, alpha, beta, Em, lam0c=None):
    m = ET.Element('material', id=str(mid), name=name, type='solid mixture')
    ET.SubElement(m, 'density').text = '1.06e-09'
    s = ET.SubElement(m, 'solid', type='neo-Hookean')
    ET.SubElement(s, 'E').text = f'{Em:.6g}'
    ET.SubElement(s, 'v').text = '0.3'
    s = ET.SubElement(m, 'solid', type='fiber-pow-linear')
    ET.SubElement(s, 'fiber', type='vector').text = '1,0,0'
    ET.SubElement(s, 'E').text = f'{E1:.6g}'
    ET.SubElement(s, 'beta').text = '2'
    ET.SubElement(s, 'lam0').text = '1.005'
    if lam0c is None:
        s = ET.SubElement(m, 'solid', type='fiber-exp-pow')
        ET.SubElement(s, 'fiber', type='vector').text = '1,0,0'
        ET.SubElement(s, 'ksi').text = f'{ksi2:.6g}'
        ET.SubElement(s, 'alpha').text = f'{alpha:.6g}'
        ET.SubElement(s, 'beta').text = f'{beta:g}'
        ET.SubElement(s, 'lam0').text = '1'
        return m
    # fiber-exp-pow-linear (FEBioMech/FEFiberPowLinear.cpp): below lam0 exactly the fitted exp-pow (its ksi from E),
    # above it linear (stress and slope continuous); E chosen so that its ksi equals the fitted ksi2
    I0 = lam0c ** 2
    E2 = ksi2 * 4 * I0 ** 1.5 * (beta - 1 + alpha * beta * (I0 - 1) ** beta) / ((I0 - 1) ** (2 - beta) * np.exp(-alpha * (I0 - 1) ** beta))
    s = ET.SubElement(m, 'solid', type='fiber-exp-pow-linear')
    ET.SubElement(s, 'fiber', type='vector').text = '1,0,0'
    ET.SubElement(s, 'E').text = f'{E2:.6g}'
    ET.SubElement(s, 'alpha').text = f'{alpha:.6g}'
    ET.SubElement(s, 'beta').text = f'{beta:g}'
    ET.SubElement(s, 'lam0').text = f'{lam0c:.6g}'
    return m


class Model109(Model7):
    def fibre_loft(self, group, m, scale=1.0, src_feb=fl.L9):
        """Swap the group's springs for its fibre lofts; matrix E = m x E1 (floor m x E1_FLOOR); fibres x scale."""
        lofts, conn_fams = GROUPS[group]
        src = Model6(src_feb)
        smesh, sdoms, scont = src.mesh, src.root.find('MeshDomains'), src.root.find('Contact')
        S, T = src.nodes(), self.nodes()
        mats = self.root.find('Material')
        md = self.root.find('MeshData')
        if md is None:
            md = ET.Element('MeshData')
            kids = list(self.root)
            self.root.insert(kids.index(self.root.find('MeshDomains')) + 1, md)
        for loft in lofts:
            blk = src.elem_blocks()[loft]
            ln = {int(v) for e in blk for v in e.text.split(',')}
            assert all(n in T and (abs(T[n] - S[n]) < 1e-12).all() for n in ln), f'{loft}: loft nodes differ'
            assert loft not in self.elem_blocks(), f'{loft} already in the target'
            _insert_like(smesh, blk, self.mesh)
            dom = copy.deepcopy(next(d for d in sdoms if d.get('name') == loft))
            E1, k2, a, b, err = law(loft, scale)
            Em = m * max(E1, E1_FLOOR)      # E1 as scaled; the floor also holds an impaired family's matrix
            mname = loft + '_fib'
            mid = max(int(x.get('id')) for x in mats) + 1
            lc = lam_cap(loft) if CAP['on'] else None
            mats.append(material_xml(mname, mid, E1, k2, a, b, Em, lc))
            if lc:
                self.log.append(f'{loft}: the stiffening fibre capped: fiber-exp-pow-linear, lam0 {lc:.3g} (the largest '
                                f'stretch of the fit range), exactly the fitted exp-pow below it, linear above')
            dom.set('mat', mname)
            doms = self.root.find('MeshDomains')
            prev = [d.get('name') for d in sdoms]
            have = [d.get('name') for d in doms]
            pos = len(doms)
            for p in reversed(prev[:prev.index(loft)]):
                if p in have:
                    pos = have.index(p) + 1
                    break
            doms.insert(pos, dom)
            # fibre directions
            eids, A, N, off = fl.fibre_dirs(fams()[loft], loft)
            order = [int(e.get('id')) for e in blk]
            assert order == eids, f'{loft}: element order differs from the survey'
            ed = ET.SubElement(md, 'ElementData', type='mat_axis', elem_set=loft)
            for i, (av, nv) in enumerate(zip(A, N), 1):
                el = ET.SubElement(ed, 'elem', lid=str(i))
                ET.SubElement(el, 'a').text = ','.join(f'{v:.6f}' for v in av)
                ET.SubElement(el, 'd').text = ','.join(f'{v:.6f}' for v in nv)
            # the loft's ties (from the source), as connectors_to_lofts
            pairs = {p.get('name'): p for p in smesh.findall('SurfacePair')}
            moved = []
            for c in scont:
                p = pairs.get(c.get('surface_pair'))
                if p is None:
                    continue
                names = [p.find('primary').text, p.find('secondary').text]
                if not any(src._surface_nodes(s) <= ln for s in names):
                    continue
                for s in smesh.findall('Surface'):
                    if s.get('name') in names:
                        _insert_like(smesh, s, self.mesh)
                _insert_like(smesh, p, self.mesh)
                _insert_like(scont, c, self.root.find('Contact'))
                moved.append(f'{c.get("name")} ({c.get("type")})')
            self.log.append(f'{LABEL}: {loft}: {len(blk)} tri3 on {len(ln)} nodes from {os.path.basename(src_feb)}, '
                            f'material {mname} = solid mixture: neo-Hookean E {Em:.4g} v 0.3 + fiber-pow-linear E1 '
                            f'{E1:.4g} (beta 2, lam0 1.005) + fiber-exp-pow ksi {k2:.4g} alpha {a:.4g} beta {b:g} '
                            f'(strip fit to its {len(fams()[loft]["rows"])} connectors, rms {100 * err:.1f} %'
                            + (f', fibres x{scale:g}' if scale != 1 else '') + f'), density 1.06e-9; fibre directions '
                            f'mat_axis (out of plane before projection: median {np.median(off):.1f}, max {off.max():.1f} deg); '
                            f'ties with it: {moved or "none"}')
        # the family's springs out
        disc = self.root.find('Discrete')
        for fam in conn_fams:
            pat = re.compile(re.escape(fam.rstrip('_')) + r'_conn(_\d+)?$')
            sets = [d for d in self.mesh.findall('DiscreteSet') if pat.match(d.get('name'))]
            assert sets, f'{fam}: no spring sets'
            gone = []
            for ds in sets:
                nm = ds.get('name')
                bb = next(x for x in disc.findall('discrete') if x.get('discrete_set') == nm)
                mat = next(x for x in disc.findall('discrete_material') if x.get('id') == bb.get('dmat'))
                assert mat.get('name') == nm + '_mat', (nm, mat.get('name'))
                disc.remove(bb)
                disc.remove(mat)
                self.mesh.remove(ds)
                gone.append(f'{nm} ({len(ds)} springs)')
            self.log.append(f'{fam}: its Abaqus connectors as springs removed: {", ".join(gone)}')
            if fam in LINE_BC:
                bc = LINE_BC[fam]
                old, new = set(self.nodesets_ids(bc)), set(src.nodesets_ids(bc))
                assert new <= old, (bc, sorted(new - old))
                cur = self._nodeset(bc)
                i = list(self.mesh).index(cur)
                self.mesh.remove(cur)
                self.mesh.insert(i, copy.deepcopy(src._nodeset(bc)))
                self.log.append(f'{bc}: back to the lofts model\'s {len(new)} nodes')
        renumber_discrete_materials(self)


TISSUE = ('_PickedSet64', '_PickedSet66', '_PickedSet346', '_PickedSet347')


HOLD = {'penalty': '10', 'maxaug': '10'}      # the edge holds' settings (batch 112: penalty 100, maxaug 0)


def _hold_edges(model, loft):
    """Hold a loft's anchor and tissue edges along their whole length (batch 110; NOT IN SOURCE, a loft convention like
    variants6.fix_anchor_edge). Walking the loft's boundary loop, every run of free nodes between two held nodes of the
    same kind is held: between two BC-fixed anchor nodes of one BC set, the run joins that set; between two chain or
    two tissue nodes, each free node follows the two by linear constraints (per dof, u - (1-s) u_A - s u_B = 0, s by arc
    length along the edge; tol 0.01, penalty 10, maxaug 10 as LA_truss_ties). Runs from an anchor to a tissue node are
    the loft's free sides and stay free. Returns (nodes added to BCs, constrained nodes)."""
    X = model.nodes()
    blocks = model.elem_blocks()
    conn = [[int(v) for v in e.text.split(',')] for e in blocks[loft]]
    owners = {}
    for nm, blk in blocks.items():
        for e in blk:
            for v in e.text.split(','):
                owners.setdefault(int(v), set()).add(nm)
    bcsets = {}          # node -> the fully fixed BC set holding it
    for bc in model.root.find('Boundary'):
        ns = bc.get('node_set')
        if not ns or bc.get('type') != 'zero displacement' or                 any((bc.findtext(d) or '0').strip() != '1' for d in ('x_dof', 'y_dof', 'z_dof')):
            continue
        try:
            ids = model.nodesets_ids(ns)
        except StopIteration:
            continue
        for n in ids:
            bcsets.setdefault(n, ns)
    chain = {int(v) for e in blocks.get('chain_mass', []) for v in e.text.split(',')}
    edges = {}
    for c in conn:
        for i in range(3):
            k = tuple(sorted((c[i], c[(i + 1) % 3])))
            edges[k] = edges.get(k, 0) + 1
    adj = {}
    for (a, b), k in edges.items():
        if k == 1:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
    start = next(iter(adj))
    loop, prev = [start], None
    while True:
        nxt = [n for n in adj[loop[-1]] if n != prev]
        if not nxt or nxt[0] == start:
            break
        prev = loop[-1]
        loop.append(nxt[0])
    assert len(loop) == len(adj), f'{loft}: boundary is not one loop ({len(loop)} of {len(adj)} nodes)'

    def kind(n):
        if n in bcsets:
            return 'bc'
        if n in chain:
            return 'chain'
        if owners[n] & set(TISSUE):
            return 'tissue'
        return None
    kinds = [kind(n) for n in loop]
    held = [i for i, k in enumerate(kinds) if k]
    added, cons = [], []
    L = len(loop)
    for j, i0 in enumerate(held):
        i1 = held[(j + 1) % len(held)]
        run = [(i0 + s) % L for s in range(1, (i1 - i0) % L)]
        side = lambda k: 'tissue' if k == 'tissue' else 'anchor'     # a BC node and a chain node: both anchors
        if not run or side(kinds[i0]) != side(kinds[i1]):
            continue
        if len(run) > 0.25 * L:
            # a run between two held nodes that goes round most of the loop is not a gap in one edge (batch 110's
            # bug: USL-L_fan, tie-attached, has only its 5 BC anchors held, and the run from the last one back to
            # the first covered its sides and tissue edge: 118 nodes were fixed)
            continue
        a, b = loop[i0], loop[i1]
        if kinds[i0] == 'bc' and kinds[i1] == 'bc' and bcsets[a] == bcsets[b]:
            added.append((bcsets[a], [loop[i] for i in run]))
            continue
        pts = [a] + [loop[i] for i in run] + [b]
        seg = np.r_[0, np.cumsum([np.linalg.norm(X[pts[q + 1]] - X[pts[q]]) for q in range(len(pts) - 1)])]
        for q, n in enumerate(pts[1:-1], 1):
            s = seg[q] / seg[-1]
            cons.append((n, a, 1 - s, b, s))
    for ns, ids in added:
        model._extend_nodeset(ns, ids)
    if cons:
        cs = model.root.find('Constraints')
        c = ET.SubElement(cs, 'constraint', name=f'{loft}_edge_holds', type='linear constraint')
        ET.SubElement(c, 'tol').text = '0.01'
        ET.SubElement(c, 'penalty').text = HOLD['penalty']
        ET.SubElement(c, 'maxaug').text = HOLD['maxaug']
        for n, a, wa, b, wb in cons:
            for dof in ('x', 'y', 'z'):
                lc = ET.SubElement(c, 'linear_constraint')
                ET.SubElement(lc, 'node', id=str(n), bc=dof).text = '1'
                ET.SubElement(lc, 'node', id=str(a), bc=dof).text = f'{-wa:.6f}'
                ET.SubElement(lc, 'node', id=str(b), bc=dof).text = f'{-wb:.6f}'
    nb = sum(len(i) for _, i in added)
    model.log.append(f'NOT IN SOURCE (loft edges held along their whole length): {loft}: boundary loop {L} nodes, held '
                     f'{len(held)} (bc {kinds.count("bc")}, chain {kinds.count("chain")}, tissue {kinds.count("tissue")}); '
                     f'{nb} free anchor-edge nodes added to their BC set ({", ".join(f"{s}: {len(i)}" for s, i in added) or "-"}); '
                     f'{len(cons)} free edge nodes held by linear constraints to their two held neighbours (constraint '
                     f'{loft}_edge_holds)' if cons else f'NOT IN SOURCE (loft edges held): {loft}: {nb} anchor-edge nodes '
                     f'added to their BC set; no constraints needed')
    return nb, len(cons)


def fibre(groups, m, scale=None, edges=False):
    scale = scale or {}

    def fn(model):
        for g in groups:
            model.fibre_loft(g, m, scale.get(g, 1.0))
        if edges:
            for g in groups:
                for loft in GROUPS[g][0]:
                    _hold_edges(model, loft)
    return fn


def paper(m, la_pct, clusl, parcus):
    """The paper's case (build_batch107.py's mapping) on the fibre lofts: the LA at la_pct % of healthy, the CL / USL
    fibre laws x clusl and the P-arcus fibre laws x parcus (the matrix scales with the fibres, above its floor)."""
    from build_batch107 import la_healthy_pct
    fn = fibre(ALL, m, {'cl': clusl, 'usl': clusl, 'parcus': parcus})

    def f(model):
        la_healthy_pct(model, la_pct)
        fn(model)
    return f


ALL = tuple(GROUPS)
BUILDS = [('L109_fibre_all_m05', BASE, fibre(ALL, 0.05)),
          ('L109_fibre_all_m02', BASE, fibre(ALL, 0.02))]
BUILDS += [(f'L109_fibre_{g}_m05', BASE, fibre((g,), 0.05)) for g in GROUPS]
# the paper's cases P1 / P2 (Luo et al. 2015) on the fibre lofts, against L107_springs_newline_rhoi0_paperP1 / P2
for _m, _tag in ((0.05, 'm05'), (0.02, 'm02')):
    BUILDS += [(f'L109_fibre_all_{_tag}_paperP1', BASE, paper(_m, 80, 0.7, 0.15)),
               (f'L109_fibre_all_{_tag}_paperP2', BASE, paper(_m, 40, 0.4, 0.15))]

if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        if os.path.exists(os.path.join(RUNS, name)):
            print(f'{name} exists; not overwriting')
            continue
        mdl = Model109(os.path.join(RUNS, base, base + '.feb'))
        fn(mdl)
        emit(name, mdl)
