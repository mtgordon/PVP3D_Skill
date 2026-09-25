"""Fourth round of single-variable variants (2026-09-23 c): the posterior-arcus -> vaginal-wall link.

Why (claude_diag/README_2026-09-23c.md): Abaqus has 26 axial CONN3D2 connectors, P-arcus-L/R-1..13, from
posterior-arcus truss nodes 1-13 to VW-PeB (behaviours LA-Y-parcus-I/II/III-x%stiff, one force-first table,
0.117 N at 4.7 mm up to 5.93 N at 47.1 mm). The FEBio model stands in for them with the lofted fans
P-arcus-L/R_fan, whose 13 vaginal-wall nodes are exactly the connector ends. The fans were attached to the
arcus through the old Posterior_Arcus_*_tube (13 shared nodes, withLA_FULL and earlier). Since the TRUSS file
replaced that tube with the spring chain, the fans' arcus edge has been free, so TRUSS_fixed2 has no
arcus-to-vaginal-wall link at all.
New operations:
  tie_parcus_fans        the 13 fan nodes per side that sat on the old tube -> the chain node at the centre
                         of their tube ring (1.00 mm away), as exact linear constraints (like the LA ties)
  add_parcus_connectors  the 26 Abaqus connectors as FEBio nonlinear springs (chain node -> PVW node),
                         force-first table swapped, elongation measure, extend constant (Abaqus default,
                         so tension-only)
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from variants import JOBS, emit, nearest_feb_dict  # noqa: E402,F401
from variants3 import Model3  # noqa: E402
from febmodel import Feb  # noqa: E402
from inpmap import part_nodes  # noqa: E402
from abq_surf import lines  # noqa: E402

WITHLA_FULL = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_FULL.feb')
SIDES = (('Left', 'L'), ('Right', 'R'))


def parcus_connectors():
    """[(name, arcus part, arcus local node, VW-PeB local node, behaviour)] and {behaviour: [(F, u)]}."""
    L = lines()
    conns = []
    for i, ln in enumerate(L):
        if ln.strip().lower().startswith('*element, type=conn3d2'):
            v = [t.strip() for t in L[i + 1].split(',')]
            m = re.search(r'elset=(P-arcus-[LR]-\d+),\s*behavior=([^,\s]+)', L[i + 2])
            if m:
                pa, na = v[1].rsplit('.', 1)
                pb, nb = v[2].rsplit('.', 1)
                conns.append((m.group(1), pa.strip('"')[:-2], int(na), int(nb), m.group(2)))
    tables = {}
    for beh in {c[4] for c in conns}:
        k = next(i for i, ln in enumerate(L) if ln.strip() == f'*Connector Behavior, name={beh}')
        assert L[k + 1].strip().lower().startswith('*connector elasticity, nonlinear, component=1'), L[k + 1]
        rows = []
        for ln in L[k + 2:]:
            if ln.strip().startswith('*'):
                break
            F, u = (float(t) for t in ln.split(',')[:2])
            rows.append((F, u))
        tables[beh] = rows
    return conns, tables


class Model4(Model3):
    def _chain_nodes(self, side):
        X = self.nodes()
        blk = next(b for b in self.mesh.findall('Nodes') if b.get('name') == f'Posterior_Arcus_{side}_centerline')
        return {int(n.get('id')): X[int(n.get('id'))] for n in blk}, X

    def tie_parcus_fans(self, ref=WITHLA_FULL):
        ref = Feb(ref)
        la = next(c for c in self.root.find('Constraints') if c.get('name') == 'LA_truss_ties')
        c = ET.SubElement(self.root.find('Constraints'), 'constraint',
                          {'name': 'Parcus_fan_ties', 'type': 'linear constraint'})
        for tag in ('tol', 'penalty', 'maxaug'):
            ET.SubElement(c, tag).text = la.find(tag).text
        pairs = []
        for side, s in SIDES:
            chain, X = self._chain_nodes(side)
            shared = set(ref.domain_nodes(f'Posterior_Arcus_{side}_tube')) & set(ref.domain_nodes(f'P-arcus-{s}_fan'))
            fan_now = {int(v) for e in self.elem_blocks()[f'P-arcus-{s}_fan'] for v in e.text.split(',')}
            assert len(shared) == 13, (side, len(shared))
            for n in sorted(shared):
                assert n in fan_now and np.allclose(X[n], ref.nodes[n]), (side, n)
                cn, d = nearest_feb_dict(chain, X[n])
                assert abs(d - 1.0) < 0.01, (side, n, cn, d)
                pairs.append((n, cn))
                for dof in 'xyz':
                    lc = ET.SubElement(c, 'linear_constraint')
                    ET.SubElement(lc, 'node', {'id': str(n), 'bc': dof}).text = '1'
                    ET.SubElement(lc, 'node', {'id': str(cn), 'bc': dof}).text = '-1'
        self.log.append(f'P-arcus fans reattached to the arcus chain: the 13 fan nodes per side shared with the old '
                        f'Posterior_Arcus_*_tube (withLA_FULL) tied to the chain node at the centre of their tube ring '
                        f'(1.00 mm), linear constraint "Parcus_fan_ties" (tol/penalty/maxaug as LA_truss_ties): {pairs}')
        return pairs

    def _fan_chain_pairs(self, ref=WITHLA_FULL):
        """{side: [(fan node, chain node)]}: the 13 fan nodes per side that sat on the old tube -> ring-centre chain node."""
        ref = Feb(ref)
        out = {}
        for side, s in SIDES:
            chain, X = self._chain_nodes(side)
            shared = set(ref.domain_nodes(f'Posterior_Arcus_{side}_tube')) & set(ref.domain_nodes(f'P-arcus-{s}_fan'))
            assert len(shared) == 13, (side, len(shared))
            out[side] = []
            for n in sorted(shared):
                assert np.allclose(X[n], ref.nodes[n]), (side, n)
                cn, d = nearest_feb_dict(chain, X[n])
                assert abs(d - 1.0) < 0.01, (side, n, cn, d)
                out[side].append((n, cn))
        return out

    def merge_parcus_fans(self):
        """Attach each P-arcus fan to the arcus chain by SHARING nodes (visible in FEBio Studio, no constraints): the
        13 fan nodes per side that sat on the old 1 mm tube are replaced in the fan's elements by the chain node at the
        centre of their ring (a 1.00 mm move), and the 8 other fan-edge nodes lying 1.00 mm off the chain (mid-segment,
        free in the original too) are moved onto the chain line so the edge stays straight. They stay unattached."""
        pairs = self._fan_chain_pairs()
        eb = self.elem_blocks()
        X = self.nodes()
        moved = []
        for side, s in SIDES:
            rep = dict(pairs[side])
            blk = eb[f'P-arcus-{s}_fan']
            n_rep = 0
            for e in blk:
                ids = [int(v) for v in e.text.split(',')]
                new = [rep.get(v, v) for v in ids]
                n_rep += sum(a != b for a, b in zip(ids, new))
                e.text = ','.join(str(v) for v in new)
            # remaining fan-edge nodes ~1 mm off the chain: project onto the chain polyline
            conns = [[int(v) for v in e.text.split(',')] for e in blk]
            ec = defaultdict(int)
            for c in conns:
                for i in range(len(c)):
                    ec[tuple(sorted((c[i], c[(i + 1) % len(c)])))] += 1
            bnd = {n for k, v in ec.items() if v == 1 for n in k}
            chain = self._chain_order(side)
            P = np.array([X[n] for n in chain])
            for n in sorted(bnd):
                if n in chain:
                    continue
                x = X[n]
                best = None
                for i in range(len(P) - 1):
                    a, b = P[i], P[i + 1]
                    tt = np.clip((x - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1)
                    q = a + tt * (b - a)
                    if best is None or np.linalg.norm(q - x) < best[0]:
                        best = (np.linalg.norm(q - x), q)
                if best[0] < 1.2:
                    self._set_node(n, best[1])
                    moved.append((n, round(best[0], 2)))
            self.log.append(f'P-arcus-{s}_fan merged onto the arcus chain: {n_rep} element-node references replaced '
                            f'(fan node -> chain node) {pairs[side]}')
        self.log.append(f'fan-edge nodes moved onto the chain line (distance mm, left unattached): {moved}')
        return pairs

    def _chain_order(self, side):
        pairs = next(ds for ds in self.mesh.findall('DiscreteSet') if ds.get('name') == f'Posterior_Arcus_{side}_springs')
        pp = [tuple(int(v) for v in e.text.split(',')) for e in pairs.findall('delem')]
        return [pp[0][0]] + [b for _, b in pp]

    def _set_node(self, nid, xyz):
        for blk in self.mesh.findall('Nodes'):
            for n in blk:
                if int(n.get('id')) == nid:
                    n.text = ','.join('%.9g' % v for v in xyz)
                    return
        raise KeyError(nid)

    def set_parcus_fan_material(self, c1=0.25908, m1=3.5, k_over_mu=50.0, name='P-arcus_fan_soft'):
        """Soft 1-term Ogden for the two P-arcus fans only, fitted so a fan pulls like its 13 Abaqus connectors
        (claude_diag/tools/fit_fan_to_connectors.py: rms 2.8 % over 2-40 mm, strip model). Density is copied from the
        fans' current material (Abaqus value in the dynamic bases). k = k_over_mu * mu0 (nearly incompressible)."""
        doms = [d for d in self.root.find('MeshDomains') if d.get('name') in ('P-arcus-L_fan', 'P-arcus-R_fan')]
        old = {d.get('mat') for d in doms}
        assert len(doms) == 2 and len(old) == 1, old
        oldname = old.pop()
        oldmat = next(m for m in self.root.find('Material') if m.get('name') == oldname)
        mats = self.root.find('Material')
        mid = max(int(m.get('id')) for m in mats) + 1
        m = ET.SubElement(mats, 'material', {'id': str(mid), 'name': name, 'type': 'Ogden'})
        ET.SubElement(m, 'density').text = oldmat.find('density').text
        mu0 = c1 / 2
        ET.SubElement(m, 'k').text = '%.6g' % (k_over_mu * mu0)
        ET.SubElement(m, 'c1').text = '%.6g' % c1
        for i in range(2, 7):
            ET.SubElement(m, f'c{i}').text = '0'
        ET.SubElement(m, 'm1').text = '%.6g' % m1
        for i in range(2, 7):
            ET.SubElement(m, f'm{i}').text = '1'
        for d in doms:
            d.set('mat', name)
        self.log.append(f'P-arcus fans only: material {oldmat.get("name")} ({oldmat.get("type")}, '
                        + ', '.join(f'{p.tag}={p.text}' for p in oldmat) + f') -> {name}: Ogden c1={c1:g}, m1={m1:g}, '
                        f'k={k_over_mu * mu0:.4g} (mu0 {mu0:.4g} MPa), density {oldmat.find("density").text}; fitted to '
                        f'the 13 LA-Y-parcus connectors per side (fit_fan_to_connectors.py)')

    def set_parcus_fan_thickness(self, t_new=0.49, name='P-arcus_fan_soft'):
        """P-arcus lofts only: shell thickness -> t_new (0.49 mm = the other lofts), with the loft material's c1 and k
        scaled by t_old/t_new so the membrane pull (the connector fit) is unchanged; the bending stiffness rises by
        (t_new/t_old)^2 against wrinkling. Call after set_parcus_fan_material (the material must be loft-only)."""
        doms = [d for d in self.root.find('MeshDomains') if d.get('name') in ('P-arcus-L_fan', 'P-arcus-R_fan')]
        assert len(doms) == 2 and all(d.get('mat') == name for d in doms), [d.get('mat') for d in doms]
        t_old = float(doms[0].find('shell_thickness').text)
        s = t_old / t_new
        for d in doms:
            d.find('shell_thickness').text = '%g' % t_new
        mat = next(m for m in self.root.find('Material') if m.get('name') == name)
        for tag in ('c1', 'k', 'density'):
            mat.find(tag).text = '%.6g' % (float(mat.find(tag).text) * s)
        self.log.append(f'P-arcus lofts: shell thickness {t_old:g} -> {t_new:g} mm (as the other lofts), c1, k and density x{s:.4g} '
                        f'(membrane pull and loft mass unchanged, bending x{(1 / s) ** 2:.3g})')

    def add_static_finish(self, t_dyn_end=1.2, static_steps=20, static_dt=0.05):
        """Two analysis steps (FEBio 4 <Step>, verified in claude_diag/mini/ms_two.feb): step 1 = the existing Control
        (dynamic ramp + settle) cut at t_dyn_end; step 2 = STATIC with loads held (load curves extend CONSTANT), same
        solver block, auto time stepper (dtmax = static_dt). Gives the exact full-load static equilibrium from the
        settled dynamic state, which a restart cannot (a restart only adds steps after the original step ends)."""
        import copy
        ctrl = self.root.find('Control')
        dt = float(ctrl.find('step_size').text)
        ctrl.find('time_steps').text = str(int(round(t_dyn_end / dt)))
        st = ET.Element('Control')
        ET.SubElement(st, 'analysis').text = 'STATIC'
        ET.SubElement(st, 'time_steps').text = str(static_steps)
        ET.SubElement(st, 'step_size').text = '%g' % static_dt
        for tag in ('plot_level', 'output_level'):
            if ctrl.find(tag) is not None:
                ET.SubElement(st, tag).text = ctrl.find(tag).text
        ts = ET.SubElement(st, 'time_stepper', {'type': 'default'})
        for tag, v in (('max_retries', 10), ('opt_iter', 15), ('dtmin', 1e-6), ('dtmax', static_dt)):
            ET.SubElement(ts, tag).text = '%g' % v
        st.append(copy.deepcopy(ctrl.find('solver')))
        kids = list(self.root)
        self.root.remove(ctrl)
        steps = ET.Element('Step')
        s1 = ET.SubElement(steps, 'step', {'id': '1', 'name': 'dynamic_ramp_settle'})
        s1.append(ctrl)
        s2 = ET.SubElement(steps, 'step', {'id': '2', 'name': 'static_finish'})
        s2.append(st)
        out = self.root.find('Output')
        self.root.insert(list(self.root).index(out), steps)
        self.log.append(f'two steps: dynamic ramp + settle cut at t = {t_dyn_end} ({ctrl.find("time_steps").text} x {dt:g}); '
                        f'then STATIC {static_steps} x {static_dt:g} with loads held (auto stepper, same solver block)')

    def truss2_to_line2(self):
        """Format only (FEBio Studio 3.2 cannot read truss2): truss2 + SolidDomain -> line2 + BeamDomain linear-truss.
        Identical solver results (claude_diag/mini/truss_mass*.feb)."""
        names = []
        for blk in self.mesh.findall('Elements'):
            if blk.get('type') == 'truss2':
                blk.set('type', 'line2')
                names.append(blk.get('name'))
        md = self.root.find('MeshDomains')
        for i, d in enumerate(list(md)):
            if d.tag == 'SolidDomain' and d.get('name') in names:
                nd = ET.Element('BeamDomain', {'name': d.get('name'), 'mat': d.get('mat'), 'type': 'linear-truss'})
                for c in d:
                    nd.append(c)
                md.remove(d)
                md.insert(i, nd)
        if names:
            self.log.append(f'format only: truss2 -> line2 + BeamDomain linear-truss for {names} (FEBio Studio)')

    def add_arcus_beam(self, E=10.0, G=4.0, A=0.01, As=5.0, I=4.0, rho=1e-9,
                       sets=('Posterior_Arcus_Left_springs', 'Posterior_Arcus_Right_springs'), clamp_ends=True):
        """NOT in the source (Abaqus T3D2 trusses have no bending stiffness): a numerical stabiliser. A
        linear-beam (Simo-Reissner) element along every posterior-arcus spring, in parallel with it, with a
        negligible axial stiffness E*A (0.1 N vs the truss's 13.9 N) so the springs keep the source's axial law.
        The line2 element is fully integrated, so a single element bends with EI + G*As*L^2/12 (shear locking,
        verified in claude_diag/mini/beam_cant.feb), and a two-spring zig-zag is resisted by G*As (needs > the
        chain compression, ~1.2 N). The chain end nodes also get zero rotation, so an end twist can't make the
        matrix singular. Verified syntax: claude_diag/mini/beam_cant.feb, beam_arc_clamped.feb."""
        pairs = []
        for ds in self.mesh.findall('DiscreteSet'):
            if ds.get('name') in sets:
                pairs += [tuple(int(v) for v in e.text.split(',')) for e in ds.findall('delem')]
        assert len(pairs) == 36, len(pairs)
        blocks = self.mesh.findall('Elements')
        eid = max(int(e.get('id')) for b in blocks for e in b) + 1
        new = ET.Element('Elements', {'type': 'line2', 'name': 'arcus_beam'})
        for i, (a, b) in enumerate(pairs):
            ET.SubElement(new, 'elem', {'id': str(eid + i)}).text = f'{a},{b}'
        self.mesh.insert(list(self.mesh).index(blocks[-1]) + 1, new)
        mats = self.root.find('Material')
        mid = max(int(m.get('id')) for m in mats) + 1
        m = ET.SubElement(mats, 'material', {'id': str(mid), 'name': 'arcus_beam_mat', 'type': 'linear-beam'})
        for tag, v in (('density', rho), ('E', E), ('G', G), ('A', A), ('A1', As), ('A2', As), ('I1', I), ('I2', I)):
            ET.SubElement(m, tag).text = '%g' % v
        ET.SubElement(self.root.find('MeshDomains'), 'BeamDomain',
                      {'name': 'arcus_beam', 'mat': 'arcus_beam_mat', 'type': 'linear-beam'})
        ends = []
        if clamp_ends:
            for name in ('BC-PosArcus_L', 'BC-PosArcus_R'):
                ends += self.nodesets_ids(name)
            ns = ET.SubElement(self.mesh, 'NodeSet', {'name': 'arcus_beam_ends'})
            ns.text = ', '.join(str(n) for n in ends)
            bc = ET.SubElement(self.root.find('Boundary'), 'bc',
                               {'name': 'arcus_beam_end_rotation', 'node_set': 'arcus_beam_ends', 'type': 'zero rotation'})
            for d in ('u_dof', 'v_dof', 'w_dof'):
                ET.SubElement(bc, d).text = '1'
        self.log.append(f'NOT IN SOURCE (numerical stabiliser): {len(pairs)} linear-beam line2 elements along the '
                        f'posterior-arcus springs {sets}, E={E:g}, G={G:g}, A={A:g} (E*A={E * A:g} N), '
                        f'A1=A2={As:g} (G*As={G * As:g} N), I1=I2={I:g} (EI={E * I:g} N mm^2), density={rho:g}; '
                        + (f'zero rotation at the chain ends {ends}' if clamp_ends else 'chain ends pinned (rotations free)'))

    def nodesets_ids(self, name):
        import re as _re
        ns = next(n for n in self.mesh.findall('NodeSet') if n.get('name') == name)
        ids = [int(v) for v in _re.split(r'[,\s]+', (ns.text or '').strip()) if v]
        return ids + [int(n.get('id')) for n in ns.findall('n')]

    def add_parcus_connectors(self, extend='constant'):
        conns, tables = parcus_connectors()
        assert len(conns) == 26, len(conns)
        eb = self.elem_blocks()
        X = self.nodes()
        hexn = {int(v) for blk in eb.values() if blk.get('type') == 'hex8' for e in blk for v in e.text.split(',')}
        PVW = {i: X[i] for i in hexn}
        chains = {side: self._chain_nodes(side)[0] for side, _ in SIDES}
        part_xyz = {}
        by_table = defaultdict(list)
        for name, part, na, nb, beh in conns:
            side = 'Left' if part.endswith('Left') else 'Right'
            if part not in part_xyz:
                part_xyz[part] = part_nodes(part)
            if 'VW-PeB' not in part_xyz:
                part_xyz['VW-PeB'] = part_nodes('VW-PeB')
            fa, da = nearest_feb_dict(chains[side], part_xyz[part][na])
            fb, db = nearest_feb_dict(PVW, part_xyz['VW-PeB'][nb])
            assert da < 1e-3 and db < 1e-3, (name, da, db)
            by_table[tuple(tables[beh])].append((name, fa, fb))
        for k, (tab, rows) in enumerate(by_table.items(), 1):
            nm = 'Parcus_conn' if len(by_table) == 1 else f'Parcus_conn_{k}'
            self.add_discrete_set(nm, [(a, b) for _, a, b in rows])
            self.add_nonlinear_spring(nm + '_mat', nm, [(u, F) for F, u in tab], extend=extend)
            self.log.append(f'added {len(rows)} Abaqus P-arcus connectors as nonlinear springs "{nm}" (elongation, '
                            f'extend={extend}; force-first table {list(tab)} swapped to (u, F)): '
                            + ', '.join(f'{n}:{a}-{b}' for n, a, b in rows))
