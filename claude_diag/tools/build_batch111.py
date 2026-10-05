"""Batch 111 (2026-09-29, the user's choice: option 3): every connector a chain of springs inside a weak sheet.
The user pictured the fibres as the source's connectors, attached at their origin point and to the tissue. In FEBio a
fibre is a direction inside each element's material (it pulls only on its element's nodes), so the fibre lofts of
batches 109-110 are held only where their edge nodes are. Here, on the springs line L87_springs_newline_rhoi0, each loft
family's connectors stay exact line elements and a sheet joins them:
  * every connector (the springs line's spring, origin node a -> tissue node b, the source's force-elongation table)
    becomes a chain of n springs along the straight line a -> b through n-1 new nodes, n = round(median length / 3 mm)
    per family (at least 2). Each segment carries the connector's table with the elongation divided by n, so a chain
    that stays straight follows the connector's law exactly (equal force, equal strain in series);
  * a sheet of tri3 shells (thickness 0.49, the lofts') between neighbouring chains of a family (ordered along the
    tissue edge, loft_fit_all.order_along_edge): a ruled surface whose every node is a chain node, so the sheet is held
    along its whole anchor and tissue edges by the connector ends; its sides are the outermost chains;
  * the sheet's material: neo-Hookean, E = m x the family's fibre E1 (build_batch109.law: the strip-model small-strain
    modulus of the connectors; floor m x 0.107 MPa), v 0.3, density 1.06e-9 (tissue); the springs are massless, as on
    the springs line.
The sheet replaces the loft surface of L9_fitall_edge (not used here): it spans only the ruled surface between the
connector lines. NOT IN SOURCE: the chains' interior nodes and the sheet (the source has single connectors).
  L111_chains_all_m05   all 14 families, sheet 5 % of E1
  L111_chains_all_m02   the same, sheet 2 %
usage: py -3.10 build_batch111.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch109 import Model109, law, BASE, RUNS, emit, E1_FLOOR  # noqa: E402
from variants2 import renumber_discrete_materials  # noqa: E402
from loft_fit_all import order_along_edge  # noqa: E402

WANT = set(sys.argv[1:])
SEG = 3.0
T_SHEET = 0.49
LABEL = 'NOT IN SOURCE (connector as a chain of springs in a sheet)'
# sheet -> (the loft whose fit gives E1, the springs line's sets, side of the P-arcus set)
SHEETS = {
    'AVW-Para-L': ('AVW-Para-L_fan', ['AVW-Para-L_conn'], None),
    'AVW-Para-R': ('AVW-Para-R_fan', ['AVW-Para-R_conn'], None),
    'CL-L': ('CL-L_fan', ['CL-L_conn_1', 'CL-L_conn_2', 'CL-L_conn_3'], None),
    'CL-R': ('CL-R_fan', ['CL-R_conn_1', 'CL-R_conn_2', 'CL-R_conn_3'], None),
    'USL-L': ('USL-L_fan', ['USL-L_conn_1', 'USL-L_conn_2'], None),
    'USL-R': ('USL-R_fan', ['USL-R_conn_1', 'USL-R_conn_2'], None),
    'PM': ('PM_fan', ['PM_conn'], None),
    'P-arcus-L': ('P-arcus-L_fan', ['Parcus_conn'], 'Posterior_Arcus_Left_springs'),
    'P-arcus-R': ('P-arcus-R_fan', ['Parcus_conn'], 'Posterior_Arcus_Right_springs'),
    'PM_PeB_Left': ('PM_PeB_Left_fan', ['PM_PeB_Left_conn'], None),
    'PM_PeB_Right': ('PM_PeB_Right_fan', ['PM_PeB_Right_conn'], None),
    'PM_avw_bottom_left': ('PM_avw_bottom_left_fan', ['PM_avw_bottom_left_conn'], None),
    'PM_avw_bottom_right': ('PM_avw_bottom_right_fan', ['PM_avw_bottom_right_conn'], None),
    'PeB-constrin': ('PeB-constrin_fan', ['PeB-constrin_conn'], None),
}


class Model111(Model109):
    def _spring_set(self, name):
        """(pairs, material element, <discrete> element, DiscreteSet element) of a spring set."""
        disc = self.root.find('Discrete')
        ds = next(d for d in self.mesh.findall('DiscreteSet') if d.get('name') == name)
        b = next(x for x in disc.findall('discrete') if x.get('discrete_set') == name)
        mat = next(x for x in disc.findall('discrete_material') if x.get('id') == b.get('dmat'))
        pairs = [tuple(int(v) for v in e.text.split(',')) for e in ds.findall('delem')]
        return pairs, mat, b, ds

    def chains(self, sheet, m, scale=1.0):
        loft, sets, side = SHEETS[sheet]
        X = self.nodes()
        side_nodes = None
        if side:
            side_nodes = {n for e in next(d for d in self.mesh.findall('DiscreteSet') if d.get('name') == side)
                          .findall('delem') for n in map(int, e.text.split(','))}
        rails = []                       # (set, a, b)
        mats = {}
        for s in sets:
            pairs, mat, b, ds = self._spring_set(s)
            mats[s] = mat
            for a, bb in pairs:
                if side_nodes is None or a in side_nodes:
                    rails.append((s, a, bb))
        assert rails, sheet
        B = np.array([X[b] for _, _, b in rails])
        order = order_along_edge(B)
        rails = [rails[i] for i in order]
        L = np.array([np.linalg.norm(X[b] - X[a]) for _, a, b in rails])
        n = max(2, int(round(np.median(L) / SEG)))
        # interior nodes of every chain
        coords, slots = [], []
        for _, a, b in rails:
            slots.append(len(coords))
            for k in range(1, n):
                coords.append(X[a] + (X[b] - X[a]) * k / n)
        ids = self.add_nodes(f'{sheet}_chain_nodes', coords)
        chain_nodes = [[a] + ids[slots[i]:slots[i] + n - 1] + [b] for i, (_, a, b) in enumerate(rails)]
        # chain springs: one set per original set, the table's elongation / n
        disc = self.root.find('Discrete')
        for s in sets:
            mem = [cn for (rs, _, _), cn in zip(rails, chain_nodes) if rs == s]
            if not mem:
                continue
            nm = f'{s}_{sheet}_chain' if side else f'{s}_chain'
            self.add_discrete_set(nm, [(cn[k], cn[k + 1]) for cn in mem for k in range(n)])
            src = mats[s]
            pts = [tuple(float(v) for v in p.text.split(',')) for p in src.find('force').find('points')]
            assert src.findtext('measure').strip() == 'elongation', (s, src.findtext('measure'))
            mid = self.add_nonlinear_spring(nm + '_mat', nm, [(u / n, F) for u, F in pts],
                                            extend=src.find('force').findtext('extend').strip())
            newm = next(x for x in disc.findall('discrete_material') if x.get('id') == str(mid))
            newm.find('scale').text = '%g' % (float(src.findtext('scale')) * scale)
        # the sheet: tri3 between neighbouring chains
        tris, skipped = [], 0
        for i in range(len(chain_nodes) - 1):
            P, Q = chain_nodes[i], chain_nodes[i + 1]
            for k in range(n):
                for tri in ((P[k], P[k + 1], Q[k + 1]), (P[k], Q[k + 1], Q[k])):
                    if len(set(tri)) < 3:
                        skipped += 1
                        continue
                    tris.append(tri)
        allX = self.nodes()
        area = [0.5 * np.linalg.norm(np.cross(allX[t[1]] - allX[t[0]], allX[t[2]] - allX[t[0]])) for t in tris]
        rung_a = [np.linalg.norm(allX[chain_nodes[i][0]] - allX[chain_nodes[i + 1][0]]) for i in range(len(rails) - 1)]
        rung_b = [np.linalg.norm(allX[chain_nodes[i][-1]] - allX[chain_nodes[i + 1][-1]]) for i in range(len(rails) - 1)]
        blocks = self.mesh.findall('Elements')
        start = max(int(e.get('id')) for blk in blocks for e in blk) + 1
        el = ET.Element('Elements', type='tri3', name=f'{sheet}_sheet')
        for q, t in enumerate(tris):
            ET.SubElement(el, 'elem', id=str(start + q)).text = ','.join(map(str, t))
        self.mesh.insert(list(self.mesh).index(blocks[-1]) + 1, el)
        E1 = law(loft, scale)[0]
        Em = m * max(E1, E1_FLOOR)
        matsec = self.root.find('Material')
        mid = max(int(x.get('id')) for x in matsec) + 1
        mat = ET.SubElement(matsec, 'material', id=str(mid), name=f'{sheet}_sheet_mat', type='neo-Hookean')
        ET.SubElement(mat, 'density').text = '1.06e-09'
        ET.SubElement(mat, 'E').text = f'{Em:.6g}'
        ET.SubElement(mat, 'v').text = '0.3'
        dom = ET.SubElement(self.root.find('MeshDomains'), 'ShellDomain', name=f'{sheet}_sheet', mat=f'{sheet}_sheet_mat')
        ET.SubElement(dom, 'shell_normal_nodal').text = '0'
        ET.SubElement(dom, 'shell_thickness').text = str(T_SHEET)
        self.log.append(f'{LABEL}: {sheet}: {len(rails)} connectors from {", ".join(sets)}'
                        + (f' ({side} side)' if side else '') + f', length {L.min():.1f}-{L.max():.1f} mm, each a chain '
                        f'of {n} springs through {n - 1} new nodes (table elongation / {n}'
                        + (f', scale x{scale:g}' if scale != 1 else '') + f'); sheet {sheet}_sheet: {len(tris)} tri3 '
                        f'({skipped} degenerate skipped), area {sum(area):.0f} mm2, smallest element {min(area):.3g} mm2, '
                        f'spacing between neighbouring chains at the anchor {min(rung_a):.2f}-{max(rung_a):.2f} mm, at '
                        f'the tissue {min(rung_b):.2f}-{max(rung_b):.2f} mm; neo-Hookean E {Em:.4g} (x{m:g} E1 '
                        f'{E1:.4g}), v 0.3, thickness {T_SHEET}')
        return sets

    def chains_all(self, m, scale=None):
        scale = scale or {}
        gone = set()
        for sheet in SHEETS:
            fam = sheet.split('-')[0] if sheet.startswith(('CL', 'USL', 'P-arcus')) else sheet
            key = {'CL': 'cl', 'USL': 'usl', 'P': 'parcus'}.get(fam, None)
            gone.update(self.chains(sheet, m, scale.get(key, 1.0)))
        # the original single springs out
        disc = self.root.find('Discrete')
        for s in sorted(gone):
            pairs, mat, b, ds = self._spring_set(s)
            disc.remove(b)
            disc.remove(mat)
            self.mesh.remove(ds)
            self.log.append(f'{s}: its {len(pairs)} single springs removed (now chains)')
        renumber_discrete_materials(self)


def build(m, scale=None, la_pct=None):
    def fn(model):
        if la_pct is not None:
            from build_batch107 import la_healthy_pct
            la_healthy_pct(model, la_pct)
        model.chains_all(m, scale)
    return fn


BUILDS = [('L111_chains_all_m05', BASE, build(0.05)),
          ('L111_chains_all_m02', BASE, build(0.02)),
          ('L111_chains_all_m05_paperP1', BASE, build(0.05, {'cl': 0.7, 'usl': 0.7, 'parcus': 0.15}, 80)),
          ('L111_chains_all_m05_paperP2', BASE, build(0.05, {'cl': 0.4, 'usl': 0.4, 'parcus': 0.15}, 40))]

if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        if os.path.exists(os.path.join(RUNS, name)):
            print(f'{name} exists; not overwriting')
            continue
        mdl = Model111(os.path.join(RUNS, base, base + '.feb'))
        fn(mdl)
        emit(name, mdl)
