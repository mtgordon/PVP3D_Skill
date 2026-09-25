"""Fifth round of single-variable variants (2026-09-24): every loft against its own Abaqus connectors.

Why (claude_diag/README_2026-09-24.md): the 12 non-P-arcus lofts all use Beam-CL-USL (isotropic elastic E 21, the
Abaqus CL/USL rigid pipe beam scaled down, density 7.8e-07). Against their connectors (tools/loft_fit_all.py) AVW-Para
is 2.4-2.5x too stiff, CL 37-58x, USL 180-350x, PM 0.6-0.8x (too soft), and five lofts (PM_PeB x2, PM_avw_bottom x2,
PeB-constrin) stand in for connectors that carry almost nothing (<= 0.003 N at 47 mm).
New operations:
  set_loft_material  the loft's own 1-term Ogden (loft_fit_all.py strip fit), k = 250 mu0, density kept
  remove_lofts       delete the loft domains (their nodes stay; unattached nodes are excluded by FEBio)
  set_loft_density   all loft materials -> one density (tissue 1.06e-9; the lofts weigh 1.65 kg against 69 g of tissue)
  set_end_time       run length only (time_steps), for screening
Added after the user's review (2026-09-24, batch 19):
  scale_surface_load   NOT IN SOURCE: Load-LA x 1/3 (the LA deforms far more than the user expects of Abaqus)
  set_tnof_penalty     the two tied-node-on-facet lofts (AVW-Para-L, USL-L) had penalty 0.0005 = loose; -> 100 N/mm
  parcus_edge_to_chain the 8 loose P-arcus loft edge nodes per side on the chain line become chain nodes
  remove_stab_springs  (variants2) with only_nodeset='BC-VW-mid': the 328 midline zero-length ground springs
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from variants import JOBS, emit  # noqa: E402,F401
from variants2 import remove_stab_springs  # noqa: E402,F401
from variants4 import Model4  # noqa: E402

NEAR_ZERO_LOFTS = ('PM_PeB_Left_fan', 'PM_PeB_Right_fan', 'PM_avw_bottom_left_fan', 'PM_avw_bottom_right_fan',
                   'PeB-constrin_fan')


class Model5(Model4):
    def _domain(self, name):
        return next(d for d in self.root.find('MeshDomains') if d.get('name') == name)

    def _material(self, name):
        return next(m for m in self.root.find('Material') if m.get('name') == name)

    def set_loft_material(self, loft, c1, m1, k_over_mu=250.0, note='', density=None):
        """The loft's own 1-term Ogden (FEBio c1 = 2 mu, m1 = alpha), k = k_over_mu * mu0, density copied from the
        loft's current material unless given (a near-zero loft needs a light one: its connectors are massless).
        Only this loft's domain changes material."""
        dom = self._domain(loft)
        oldmat = self._material(dom.get('mat'))
        mats = self.root.find('Material')
        name = loft + '_fit'
        assert not any(m.get('name') == name for m in mats), name
        mid = max(int(m.get('id')) for m in mats) + 1
        m = ET.SubElement(mats, 'material', {'id': str(mid), 'name': name, 'type': 'Ogden'})
        ET.SubElement(m, 'density').text = oldmat.find('density').text if density is None else '%g' % density
        mu0 = c1 / 2
        ET.SubElement(m, 'k').text = '%.6g' % (k_over_mu * mu0)
        ET.SubElement(m, 'c1').text = '%.6g' % c1
        for i in range(2, 7):
            ET.SubElement(m, f'c{i}').text = '0'
        ET.SubElement(m, 'm1').text = '%.6g' % m1
        for i in range(2, 7):
            ET.SubElement(m, f'm{i}').text = '1'
        dom.set('mat', name)
        self.log.append(f'{loft} only: material {oldmat.get("name")} ({oldmat.get("type")}, '
                        + ', '.join(f'{p.tag}={p.text}' for p in oldmat) + f') -> {name}: Ogden c1={c1:.6g}, m1={m1:g}, '
                        f'k={k_over_mu * mu0:.6g} ({k_over_mu:g} mu0), density {m.find("density").text}; {note}')

    def set_solver(self, **kw):
        """Solver-control tags by name: time_stepper (opt_iter, cutback, dtmax, max_retries), solver (lstol, max_refs,
        ...), qn_method (max_ups; type= via qn_type)."""
        ctrl = self.root.find('Control')
        ts, sv = ctrl.find('time_stepper'), ctrl.find('solver')
        qn = sv.find('qn_method')
        for key, val in kw.items():
            if key == 'qn_type':
                old = qn.get('type')
                qn.set('type', val)
            else:
                el = next((p.find(key) for p in (ts, sv, qn) if p is not None and p.find(key) is not None), None)
                assert el is not None, key
                old = el.text
                el.text = str(val)
            self.log.append(f'solver: {key} {old} -> {val}')

    def remove_lofts(self, names=NEAR_ZERO_LOFTS):
        for n in names:
            self._domain(n)       # must exist
        self.remove_domains(list(names))

    def set_loft_density(self, rho=1.06e-9):
        md = self.root.find('MeshDomains')
        mats = {d.get('mat') for d in md if d.get('name', '').endswith('_fan')}
        others = [d.get('name') for d in md if d.get('mat') in mats and not d.get('name', '').endswith('_fan')]
        assert not others, others
        for mn in sorted(mats):
            m = self._material(mn)
            old = m.find('density').text
            m.find('density').text = '%g' % rho
            self.log.append(f'loft material {mn}: density {old} -> {rho:g} (tissue)')

    def scale_surface_load(self, name='Load-LA', factor=1 / 3):
        """NOT IN SOURCE (user's choice, 2026-09-24): scale one pressure load (Abaqus: 0.014 MPa on all five)."""
        sl = next(s for s in self.root.find('Loads') if s.get('name') == name)
        p = sl.find('pressure')
        old = float(p.text)
        p.text = '%.6g' % (old * factor)
        self.log.append(f'NOT IN SOURCE: {name} pressure {old:g} -> {old * factor:.6g} MPa (x{factor:.4g}); '
                        f'the other pressure loads unchanged')

    def set_tnof_penalty(self, penalty=100.0):
        """Every tied-node-on-facet contact (AVW-Para-L and USL-L lofts) -> penalty [N/mm per node]. 0.0005 left the
        ties loose (mini/tnof_pen_*: gap = the full pull; 10 -> 0.02 mm, 100 -> 0.003 mm at ~0.3 N per node)."""
        n = 0
        for c in self.root.find('Contact'):
            if c.get('type') == 'tied-node-on-facet':
                old = c.find('penalty').text
                c.find('penalty').text = '%g' % penalty
                self.log.append(f'tied-node-on-facet {c.get("name")}: penalty {old} -> {penalty:g}')
                n += 1
        assert n, 'no tied-node-on-facet contacts'

    def parcus_edge_to_chain(self, tol=0.05):
        """The P-arcus loft's edge nodes that lie on the arcus chain line but are not chain nodes (8 per side after
        merge_parcus_fans) become chain nodes: the chain spring they lie on is split through them, in order. The
        springs are strain-based, so each piece keeps the same force-strain law; the chain mass elements stay as
        they are (the loft already gives these nodes mass)."""
        X = self.nodes()
        els = self.elem_blocks()
        for side, s in (('Left', 'L'), ('Right', 'R')):
            loft = els[f'P-arcus-{s}_fan']
            conn = [[int(v) for v in e.text.split(',')] for e in loft]
            cnt = {}
            for c in conn:
                for k in range(len(c)):
                    e = tuple(sorted((c[k], c[(k + 1) % len(c)])))
                    cnt[e] = cnt.get(e, 0) + 1
            bnd = {n for e, k in cnt.items() if k == 1 for n in e}
            ds = next(d for d in self.mesh.findall('DiscreteSet') if d.get('name') == f'Posterior_Arcus_{side}_springs')
            springs = [tuple(int(v) for v in d.text.split(',')) for d in ds.findall('delem')]
            chain_nodes = {n for sp in springs for n in sp}
            on = {}
            for n in sorted(bnd - chain_nodes):
                best = None
                for i, (a, b) in enumerate(springs):
                    A, B = X[a], X[b]
                    d = B - A
                    t = float((X[n] - A) @ d / (d @ d))
                    dist = float(np.linalg.norm(A + min(max(t, 0), 1) * d - X[n]))
                    if 0 < t < 1 and dist < tol and (best is None or dist < best[0]):
                        best = (dist, i, t)
                if best:
                    on.setdefault(best[1], []).append((best[2], n, best[0]))
            new = []
            for i, (a, b) in enumerate(springs):
                pts = [a] + [n for t, n, _ in sorted(on.get(i, []))] + [b]
                new.extend(zip(pts[:-1], pts[1:]))
            for d in ds.findall('delem'):
                ds.remove(d)
            ds.text = '\n\t\t\t'
            for a, b in new:
                d = ET.SubElement(ds, 'delem')
                d.text = f'{a},{b}'
                d.tail = '\n\t\t\t'
            added = sorted((n, round(dist, 3)) for v in on.values() for t, n, dist in v)
            self.log.append(f'P-arcus-{s} loft edge nodes on the chain line made chain nodes ({len(added)}; node, distance '
                            f'mm): {added}; Posterior_Arcus_{side}_springs {len(springs)} -> {len(new)} springs')

    def set_end_time(self, t_end):
        ctrl = self.root.find('Control')
        dt = float(ctrl.find('step_size').text)
        n = int(round(t_end / dt))
        old = ctrl.find('time_steps').text
        ctrl.find('time_steps').text = str(n)
        self.log.append(f'run length only: time_steps {old} -> {n} (ends at t = {n * dt:g}); the model is unchanged')
