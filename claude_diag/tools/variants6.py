"""Sixth round (2026-09-24, task 3): the in-situ check of the fitted lofts against the connectors they stand in for.

Why (SKILL item 12, claude_diag/README_2026-09-24.md): a strip fit is good to ~+-50 %; the exact reference is a run with
the Abaqus CONN3D2 connectors themselves as FEBio springs, whose connector-end elongations the fitted lofts should
reproduce (done for P-arcus in batch 17: LPC_DM1_pconn).
New operation:
  lofts_to_connectors  for each named connector family (loft_survey.py names): its connectors as FEBio nonlinear
                       springs (elongation; the force-first table swapped to (u, F); extend constant = the Abaqus
                       default EXTRAPOLATION=CONSTANT, so no force in compression), one DiscreteSet per family and
                       distinct table; the family's loft domain removed, with every contact whose surface lies on it
                       (the tied-node-on-facet ties of AVW-Para-L and USL-L). Ends by exact coordinate: the tissue end
                       on the VW-PeB solid node; the anchor end on the BC-fixed Abaqus RP node (AVW-Para, PM, PM_PeB,
                       PM_avw_bottom, PeB-constrin), the arcus chain node (P-arcus), or the CL/USL line node. Abaqus
                       holds each CL/USL line fixed (*Rigid Body on RP 1-4, which BC-CL/USL-Left/Right fix in all 6
                       DOFs); in FEBio those line nodes are orphans, so they are added to BC-CL/USL-Left/Right.
  fix_anchor_edge      (user, 2026-09-24) a loft's whole anchor edge under its BC: the free loft nodes between the
                       BC-fixed nodes along the loft boundary join the BC's node set. PM_PeB: Abaqus pins only the 8
                       connector origins per side; the loft edge between them had 7 (left) / 9 (right) free midpoint
                       nodes, which barely move with the E 21 loft (<= 0.15 mm) but swing up to 3.1-3.6 mm with the
                       weak fitted loft (L6_pmpeb). A loft convention, NOT IN SOURCE (the source has no loft).
"""
import os
import sys
import xml.etree.ElementTree as ET
from collections import OrderedDict

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from variants import JOBS, emit  # noqa: E402,F401
from variants5 import Model5, NEAR_ZERO_LOFTS  # noqa: E402,F401
from loft_survey import survey, FAMILY_LOFT  # noqa: E402

LINE_BC = {'CL-L': 'BC-CL-Left', 'CL-R': 'BC-CL-Right', 'USL-L': 'BC-USL-Left', 'USL-R': 'BC-USL-Right'}
# every non-P-arcus loft family (PM-middle is already the FEBio spring PM-LA-x%stiff_mat)
LOFT_FAMILIES = ('AVW-Para-L', 'AVW-Para-R', 'CL-L', 'CL-R', 'USL-L', 'USL-R', 'PM', 'PM_PeB_Left_',
                 'PM_PeB_Right_', 'PM_avw_bottom_left_', 'PM_avw_bottom_right_', 'PeB-constrin')


class Model6(Model5):
    def _boundary_loop(self, loft):
        """The loft's boundary nodes in loop order (one loop expected)."""
        conns = [[int(v) for v in e.text.split(',')] for e in self.elem_blocks()[loft]]
        cnt = {}
        for c in conns:
            for k in range(len(c)):
                e = tuple(sorted((c[k], c[(k + 1) % len(c)])))
                cnt[e] = cnt.get(e, 0) + 1
        adj = {}
        for (a, b), k in cnt.items():
            if k == 1:
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)
        assert all(len(v) == 2 for v in adj.values()), f'{loft}: boundary is not a simple loop'
        s = min(adj)
        loop, prev, cur = [s], None, s
        while True:
            nxt = [n for n in adj[cur] if n != prev and n != s]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            loop.append(cur)
        assert len(loop) == len(adj), f'{loft}: more than one boundary loop'
        return loop

    def fix_anchor_edge(self, loft, bc):
        """Every loft node on the anchor edge -> BC node set `bc`: the anchor edge is the arc of the boundary loop that
        spans all the loft's nodes in `bc` (the complement of the largest gap between them); its nodes that are not
        in `bc` must belong to this loft only (no tissue or other-domain nodes), and are added."""
        loop = self._boundary_loop(loft)
        fixed = set(self.nodesets_ids(bc))
        I = [i for i, n in enumerate(loop) if n in fixed]
        assert len(I) >= 2, (loft, bc, len(I))
        gaps = [((I[(k + 1) % len(I)] - I[k]) % len(loop), k) for k in range(len(I))]
        big = max(gaps)[1]
        start = I[(big + 1) % len(I)]
        arc = [loop[(start + j) % len(loop)] for j in range((I[big] - start) % len(loop) + 1)]
        others = {}
        for name, blk in self.elem_blocks().items():
            if name != loft:
                for e in blk:
                    for v in e.text.split(','):
                        others.setdefault(int(v), name)
        add = [n for n in arc if n not in fixed]
        bad = [(n, others[n]) for n in add if n in others]
        assert not bad, f'{loft}: anchor-edge nodes also in other domains {bad}'
        X = self.nodes()
        # distance of each added node to the straight segment between its fixed neighbours on the arc
        dist = []
        for n in add:
            j = arc.index(n)
            a = next(arc[i] for i in range(j, -1, -1) if arc[i] in fixed)
            b = next(arc[i] for i in range(j, len(arc)) if arc[i] in fixed)
            A, B, P = X[a], X[b], X[n]
            t = float(np.clip((P - A) @ (B - A) / ((B - A) @ (B - A)), 0, 1))
            dist.append(float(np.linalg.norm(A + t * (B - A) - P)))
        new = self._extend_nodeset(bc, add)
        self.log.append(f'NOT IN SOURCE (loft convention, user): {loft} whole anchor edge under {bc}: + {len(new)} free edge '
                        f'nodes between the {len(I)} fixed ones {new} (off the straight segments between them: max '
                        f'{max(dist) if dist else 0:.3g} mm)')
        return new

    def set_loft_ogden(self, loft, terms, k_over_mu=250.0, density=None, suffix='_sec', note=''):
        """The loft's own N-term Ogden, terms = [(c_i, m_i)] (FEBio c_i = 2 mu_i, m_i = alpha_i; mu0 = sum c_i / 2),
        k = k_over_mu * mu0; density copied from the loft's current material unless given."""
        dom = self._domain(loft)
        oldmat = self._material(dom.get('mat'))
        mats = self.root.find('Material')
        name = loft + suffix
        assert not any(m.get('name') == name for m in mats), name
        assert len(terms) <= 6 and all(c * m > 0 for c, m in terms), terms
        mid = max(int(m.get('id')) for m in mats) + 1
        m = ET.SubElement(mats, 'material', {'id': str(mid), 'name': name, 'type': 'Ogden'})
        ET.SubElement(m, 'density').text = oldmat.find('density').text if density is None else '%g' % density
        mu0 = sum(c for c, _ in terms) / 2
        ET.SubElement(m, 'k').text = '%.6g' % (k_over_mu * mu0)
        cs = [c for c, _ in terms] + [0] * (6 - len(terms))
        ms = [mm for _, mm in terms] + [1] * (6 - len(terms))
        for i, c in enumerate(cs, 1):
            ET.SubElement(m, f'c{i}').text = '%.6g' % c
        for i, mm in enumerate(ms, 1):
            ET.SubElement(m, f'm{i}').text = '%.6g' % mm
        dom.set('mat', name)
        self.log.append(f'{loft} only: material {oldmat.get("name")} -> {name}: Ogden '
                        + ' + '.join(f'(c {c:.4g}, m {mm:g})' for c, mm in terms)
                        + f', k={k_over_mu * mu0:.4g} ({k_over_mu:g} mu0, mu0 {mu0:.4g}), density {m.find("density").text}; {note}')

    def set_lofts_density(self, lofts, rho=1.06e-9):
        """Named lofts' own materials -> density rho (each loft must have a material of its own; Model5.set_loft_density
        sets every loft material at once)."""
        for loft in lofts:
            mn = self._domain(loft).get('mat')
            users = [d.get('name') for d in self.root.find('MeshDomains') if d.get('mat') == mn]
            assert users == [loft], (loft, users)
            mat = self._material(mn)
            old = mat.find('density').text
            mat.find('density').text = '%g' % rho
            self.log.append(f'{loft}: material {mn} density {old} -> {rho:g} (tissue)')

    def _surface_nodes(self, name):
        s = next(s for s in self.mesh.findall('Surface') if s.get('name') == name)
        return {int(v) for f in s for v in f.text.split(',')}

    def remove_contacts_on(self, nodes):
        """Remove every contact (with its SurfacePair and both Surfaces) one of whose surfaces lies entirely on the
        given nodes."""
        pairs = {p.get('name'): p for p in self.mesh.findall('SurfacePair')}
        gone = []
        for c in list(self.root.find('Contact')):
            p = pairs.get(c.get('surface_pair'))
            if p is None:
                continue
            surfs = [p.find('primary').text, p.find('secondary').text]
            if any(self._surface_nodes(s) <= nodes for s in surfs):
                self.root.find('Contact').remove(c)
                self.mesh.remove(p)
                for s in list(self.mesh.findall('Surface')):
                    if s.get('name') in surfs:
                        self.mesh.remove(s)
                gone.append(c.get('name'))
        return gone

    def lofts_to_connectors(self, families=LOFT_FAMILIES):
        f, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(self.src)
        dm = self.root.find('Discrete').findall('discrete_material')
        assert [int(m.get('id')) for m in dm] == list(range(1, len(dm) + 1)), 'dmat ids must be list positions'
        ids = np.array(sorted(f.nodes))
        P = np.array([f.nodes[i] for i in ids])
        solid = {n for et, d in f.elem_blocks.values() if et == 'hex8' for c in d.values() for n in c}
        chain = {n for nm, prs in f.discsets.items() if nm.startswith('Posterior_Arcus_') for p in prs for n in p}
        lofts = sorted({FAMILY_LOFT[fam] for fam in families})
        loft_nodes = {n for lf in lofts for c in f.elem_blocks[lf][1].values() for n in c}

        def at(xyz):
            return [int(n) for n in ids[np.linalg.norm(P - xyz, axis=1) < 1e-6]]

        fixed = {}
        for fam in families:
            loft, rows = res[fam]
            assert loft in lofts, (fam, loft)
            groups = OrderedDict()
            for r in rows:
                ca = at(r['xa'])
                pick = ([n for n in ca if n in bcs] or [n for n in ca if n in chain]
                        or [n for n in ca if not (elem_of[n] - set(lofts))])
                assert pick, (fam, r['elset'], ca, [sorted(elem_of[n]) for n in ca])
                a = min(pick)
                if a not in bcs and a not in chain:
                    assert fam in LINE_BC, (fam, r['elset'], a, 'anchor neither fixed nor on the chain')
                    fixed.setdefault(LINE_BC[fam], set()).add(a)
                cb = [n for n in at(r['xb']) if n in solid]
                assert cb, (fam, r['elset'], 'no VW-PeB solid node at the tissue end')
                tab = tuple(map(tuple, behav[r['behavior']][0]['table']))
                groups.setdefault(tab, []).append((r['elset'], a, min(cb)))
            for k, (tab, members) in enumerate(groups.items(), 1):
                nm = fam.rstrip('_') + '_conn' + ('' if len(groups) == 1 else f'_{k}')
                self.add_discrete_set(nm, [(a, b) for _, a, b in members])
                self.add_nonlinear_spring(nm + '_mat', nm, [(u, F) for F, u in tab], extend='constant')
                self.log.append(f'{fam}: {len(members)} Abaqus connectors as nonlinear springs "{nm}" (elongation, '
                                f'extend constant; force-first table {list(tab)} swapped to (u, F)): '
                                + ', '.join(f'{e}:{a}-{b}' for e, a, b in members))
        for bc, nodes in fixed.items():
            new = self._extend_nodeset(bc, sorted(nodes))
            self.log.append(f'{bc}: + {len(new)} CL/USL line nodes (connector anchors; Abaqus: rigid body on a fully '
                            f'fixed RP) {new}')
        gone = self.remove_contacts_on(loft_nodes)
        self.remove_domains(lofts)
        self.log.append(f'contacts on the removed lofts removed (with their surface pairs and surfaces): {gone}')
