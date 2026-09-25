"""Third round of single-variable variants (2026-09-23 b): LA material fidelity + Load-LA strategies.

Base model: PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed2.feb (= run ALL2_faithful_fans).
Why (claude_diag/la_stiffness/): the Abaqus LA material LA_Yamada50% is Marlow, whose deviatoric energy depends
on I1 only. The FEBio LA is a 1-term Ogden fitted to the same uniaxial data: equal in uniaxial tension, but
2-3x softer in equibiaxial stretch at 20-40 % strain. A 2-term Yeoh (also I1-only) fitted to the same data
reproduces Marlow within ~10 % in uniaxial, planar and equibiaxial stretch. Separately, k from mu0 (0.33)
lets the FEBio LA lose volume stiffness as it stiffens; Abaqus keeps nu = 0.47 at every strain, which needs
k ~ 0.5-1.3 MPa over 20-50 % strain.
New operations:
  set_la_yeoh     LA material -> uncoupled Yeoh (c1, c2, k)
  set_la_k        LA bulk modulus
  delay_la_load   Load-LA on its own smooth-step curve from t0 to t1 (other loads unchanged)
"""
import os
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(__file__))
from variants import Model, JOBS, emit  # noqa: E402,F401
from variants2 import Model2  # noqa: E402

TRUSS_FIXED2 = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed2.feb')
YEOH_C = (0.0100393, 0.0172014)   # fit to Abaqus LA_Yamada50% uniaxial data (la_stiffness/la_yeoh_fit.py)
K_MU0 = 0.32795                   # k = 2 mu0 (1+nu) / (3 (1-2nu)), mu0 = 2 c1, nu = 0.47


CHAIN_SETS = ('ATLA_Left_springs', 'ATLA_Right_springs', 'Posterior_Arcus_Left_springs',
              'Posterior_Arcus_Right_springs')


class Model3(Model2):
    ABQ_DENSITY = {**Model.ABQ_DENSITY, 'LA_Yamada50pct_Yeoh': 1.06e-09, 'chain_mass_mat': 0.00011}

    def add_chain_mass(self, sets=CHAIN_SETS, rho=0.00011, area=1.0, E=1e-6):
        """Abaqus ATLA / Posterior_Arcus trusses: T3D2, area 1, *Density 0.00011 (~1e5 x tissue). The FEBio
        chains are massless springs, so for a DYNAMIC run add a truss element (line2, BeamDomain linear-truss) along every chain spring with
        that density and a negligible modulus (mass only). Verified syntax: claude_diag/mini/truss_mass.feb
        (consistent mass rho*A*L/3 at the free node)."""
        pairs = []
        for ds in self.mesh.findall('DiscreteSet'):
            if ds.get('name') in sets:
                pairs += [e.text.strip() for e in ds.findall('delem')]
        blocks = self.mesh.findall('Elements')
        eid = max(int(e.get('id')) for b in blocks for e in b) + 1
        # line2 + BeamDomain linear-truss (documented FEBio 4 form; FEBio Studio 3.2 can't read truss2, and the
        # solver gives identical results: claude_diag/mini/truss_mass.feb vs truss_mass_line2.feb)
        new = ET.Element('Elements', {'type': 'line2', 'name': 'chain_mass'})
        for i, p in enumerate(pairs):
            ET.SubElement(new, 'elem', {'id': str(eid + i)}).text = p
        self.mesh.insert(list(self.mesh).index(blocks[-1]) + 1, new)
        mats = self.root.find('Material')
        mid = max(int(m.get('id')) for m in mats) + 1
        m = ET.SubElement(mats, 'material', {'id': str(mid), 'name': 'chain_mass_mat', 'type': 'linear truss'})
        ET.SubElement(m, 'density').text = '%g' % rho
        ET.SubElement(m, 'E').text = '%g' % E
        ET.SubElement(m, 'v').text = '0.3'
        dom = ET.SubElement(self.root.find('MeshDomains'), 'BeamDomain',
                            {'name': 'chain_mass', 'mat': 'chain_mass_mat', 'type': 'linear-truss'})
        ET.SubElement(dom, 'cross_sectional_area').text = '%g' % area
        self.log.append(f'added {len(pairs)} mass-only truss elements (line2, BeamDomain linear-truss) along the chain springs {sets} '
                        f'(linear truss, rho={rho:g}, E={E:g}, area={area:g}; Abaqus truss density)')

    def _la_material(self):
        doms = [d for d in self.root.find('MeshDomains') if d.get('name', '').startswith('LA_')]
        names = {d.get('mat') for d in doms}
        assert len(names) == 1, names
        name = names.pop()
        return next(m for m in self.root.find('Material') if m.get('name') == name), doms

    def set_la_yeoh(self, c=YEOH_C, k=K_MU0, name='LA_Yamada50pct_Yeoh'):
        mat, doms = self._la_material()
        old = f"{mat.get('type')} " + ', '.join(f'{p.tag}={p.text}' for p in mat)
        dens = mat.find('density').text
        for p in list(mat):
            mat.remove(p)
        mat.set('type', 'Yeoh')
        mat.set('name', name)
        ET.SubElement(mat, 'density').text = dens
        for i, ci in enumerate(c, 1):
            ET.SubElement(mat, f'c{i}').text = '%.7g' % ci
        ET.SubElement(mat, 'k').text = '%.6g' % k
        for d in doms:
            d.set('mat', name)
        self.log.append(f'LA material -> uncoupled Yeoh (I1-only, like Abaqus Marlow LA_Yamada50%): '
                        f'c={c}, k={k:g} (was {old})')

    def set_la_k(self, k):
        mat, _ = self._la_material()
        old = mat.find('k').text
        mat.find('k').text = '%.6g' % k
        self.log.append(f'LA bulk modulus k {old} -> {k:g}')

    def scale_la(self, factor):
        """NOT faithful: LA stiffer than Abaqus (every c_i and k times factor)."""
        mat, _ = self._la_material()
        for p in mat:
            if p.tag in ('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'k'):
                p.text = '%.6g' % (float(p.text) * factor)
        self.log.append(f'NOT FAITHFUL: LA material c_i and k x{factor:g} (stiffer than Abaqus)')

    def add_settle_phase(self, t_end=2.0, C=20.0, t_on=1.0, t_full=1.05):
        """Dynamic relaxation after the ramp: run on to t_end with every load held (load curves extend
        CONSTANT) and mass-proportional damping C [1/s] switched on from t_on (0 before, so the ramp is
        unchanged). The settled state is the static equilibrium at full load, whatever the masses.
        Verified syntax: claude_diag/mini/truss_massdamp.feb."""
        ctrl = self.root.find('Control')
        dt = float(ctrl.find('step_size').text)
        ctrl.find('time_steps').text = str(int(round(t_end / dt)))
        ld = self.root.find('LoadData')
        nid = max(int(lc.get('id')) for lc in ld) + 1
        lc = ET.SubElement(ld, 'load_controller', {'id': str(nid), 'name': 'settle_on', 'type': 'loadcurve'})
        ET.SubElement(lc, 'interpolate').text = 'SMOOTH STEP'
        ET.SubElement(lc, 'extend').text = 'CONSTANT'
        pts = ET.SubElement(lc, 'points')
        ET.SubElement(pts, 'pt').text = f'{t_on},0'
        ET.SubElement(pts, 'pt').text = f'{t_full},1'
        bl = ET.SubElement(self.root.find('Loads'), 'body_load', {'name': 'settle_damping', 'type': 'mass damping'})
        ET.SubElement(bl, 'C', {'lc': str(nid)}).text = '%g' % C
        self.log.append(f'settle phase: run to t={t_end} with loads held, mass damping C={C:g}/s ramped in '
                        f'over t={t_on}-{t_full} (lc {nid})')

    def delay_la_load(self, t0=0.5, t1=1.0):
        ld = self.root.find('LoadData')
        nid = max(int(lc.get('id')) for lc in ld) + 1
        lc = ET.SubElement(ld, 'load_controller', {'id': str(nid), 'name': 'LA_late', 'type': 'loadcurve'})
        ET.SubElement(lc, 'interpolate').text = 'SMOOTH STEP'
        ET.SubElement(lc, 'extend').text = 'CONSTANT'
        pts = ET.SubElement(lc, 'points')
        ET.SubElement(pts, 'pt').text = f'{t0},0'
        ET.SubElement(pts, 'pt').text = f'{t1},1'
        n = 0
        for sl in self.root.find('Loads'):
            if sl.get('name') == 'Load-LA':
                sl.find('pressure').set('lc', str(nid))
                n += 1
        assert n == 1, n
        self.log.append(f'Load-LA on its own smooth-step curve (lc {nid}): 0 until t={t0}, full at t={t1}')
