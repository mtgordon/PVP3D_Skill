"""Second round of single-variable variants (2026-09-23), built on variants.Model.

Base model: PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed.feb (= run T3f_k_snn0).
New operations:
  add_pvw_la_contact  Abaqus Int-PVW-LA (frictionless, PVW back vs LA SPOS) as sliding-elastic
  add_la_pressure     Abaqus Load-LA (0.014 MPa on the SPOS face of all 1170 LA elements, Amp-1)
  add_vw_mid_xsymm    Abaqus BC-VW-mid (XSYMM on _PickedSet1359 = FEBio NodeSet BC-VW-mid)
  set_connector_extend  FEBio spring curve extension for the Abaqus-connector springs
  add_peb_bottom_to_avw_contact  the 208 PeB-bottom facets Abaqus has in Int-AVW-PVW-PEB-inner
  set_diagnostic_output  per-iteration plotting of displacement + relative volume only
"""
import os
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(__file__))
from variants import Model, JOBS, emit  # noqa: E402,F401

BEST = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed.feb')
LA_DOMAINS = ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri')
CONNECTOR_SPRINGS = ('LA_sphincter_side_conn_mat', 'LA_sphincter_post_conn_mat', 'PM-LA-x%stiff_mat')


class Model2(Model):
    # ---------- helpers ----------
    def _add_surface(self, name, facets):
        """facets: list of node-id lists (3 or 4 nodes), appended at the end of <Mesh>."""
        s = ET.Element('Surface', {'name': name})
        for i, c in enumerate(facets, 1):
            f = ET.SubElement(s, 'quad4' if len(c) == 4 else 'tri3', {'id': str(i)})
            f.text = ','.join(str(n) for n in c)
        self.mesh.append(s)          # end of <Mesh>: after every <Nodes> block it references (gotcha 1)
        return s

    def _add_surface_pair(self, name, primary, secondary):
        sp = ET.Element('SurfacePair', {'name': name})
        ET.SubElement(sp, 'primary').text = primary
        ET.SubElement(sp, 'secondary').text = secondary
        self.mesh.append(sp)

    def la_facets(self):
        eb = self.elem_blocks()
        return [[int(v) for v in e.text.split(',')] for n in LA_DOMAINS for e in eb[n]]

    def surface_facets(self, name):
        for s in self.mesh.findall('Surface'):
            if s.get('name') == name:
                return [[int(v) for v in f.text.split(',')] for f in s]
        raise KeyError(name)

    # ---------- operations ----------
    def add_pvw_la_contact(self, offset=2.0, penalty=5, search_radius=3.0, two_pass=1, auto_penalty=1,
                           laugon='PENALTY'):
        """Abaqus *Contact Pair Int-PVW-LA (Frictionless, KINEMATIC): Surf-VW-PVW-back vs _PickedSurf479
        (SPOS of all LA elements). FEBio LA shell nodes are the shell TOP face (+n side, the side the PVW
        is on); the Abaqus SPOS face is 2 mm (t/2) above the nodes, hence offset=2. Parameters otherwise
        copy the model's validated SlidingElastic1 (AVW-PVW) contact."""
        pv = self.surface_facets('Load-PVW_surf')          # == Abaqus Surf-VW-PVW-back (676/676 by coordinate)
        self._add_surface('PVW_LA_primary', pv)
        self._add_surface('PVW_LA_secondary', self.la_facets())
        self._add_surface_pair('PVW_LA', 'PVW_LA_primary', 'PVW_LA_secondary')
        con = self.root.find('Contact')
        c = ET.SubElement(con, 'contact', {'name': 'PVW_LA', 'surface_pair': 'PVW_LA', 'type': 'sliding-elastic'})
        for k, v in (('laugon', laugon), ('tolerance', 0.1), ('gaptol', 0), ('penalty', penalty),
                     ('auto_penalty', auto_penalty), ('update_penalty', 0), ('two_pass', two_pass),
                     ('knmult', 0), ('search_tol', 0.01), ('symmetric_stiffness', 0),
                     ('search_radius', search_radius), ('seg_up', 0), ('tension', 0), ('minaug', 0),
                     ('maxaug', 15), ('node_reloc', 0), ('fric_coeff', 0), ('smooth_aug', 0),
                     ('flip_primary', 0), ('flip_secondary', 0), ('shell_bottom_primary', 0),
                     ('shell_bottom_secondary', 0), ('offset', offset)):
            ET.SubElement(c, k).text = str(v)
        self.log.append(f'added Abaqus Int-PVW-LA contact: sliding-elastic PVW back (676, = Load-PVW_surf) vs '
                        f'LA top face (all {len(self.la_facets())} LA elements), offset={offset}, penalty={penalty} '
                        f'(auto={auto_penalty}), two_pass={two_pass}, search_radius={search_radius}, {laugon}')

    def add_la_pressure(self, p=0.014, lc='1'):
        """Abaqus Load-LA: *Dsload _PickedSurf484 (SPOS of all 1170 LA elements), P=0.014, Amp-1.
        FEBio pressure on the element-ordered LA facets acts along -n, like Abaqus P on SPOS."""
        self._add_surface('Load-LA_surf', self.la_facets())
        loads = self.root.find('Loads')
        sl = ET.SubElement(loads, 'surface_load', {'name': 'Load-LA', 'surface': 'Load-LA_surf', 'type': 'pressure'})
        ET.SubElement(sl, 'pressure', {'lc': lc}).text = str(p)
        ET.SubElement(sl, 'symmetric_stiffness').text = '1'
        ET.SubElement(sl, 'linear').text = '0'
        ET.SubElement(sl, 'shell_bottom').text = '0'
        self.log.append(f'added Abaqus Load-LA: pressure {p} (lc {lc}) on the top face of all LA elements')

    def add_vw_mid_xsymm(self):
        bnd = self.root.find('Boundary')
        bc = ET.SubElement(bnd, 'bc', {'name': 'BC-VW-mid', 'node_set': 'BC-VW-mid', 'type': 'zero displacement'})
        ET.SubElement(bc, 'x_dof').text = '1'
        ET.SubElement(bc, 'y_dof').text = '0'
        ET.SubElement(bc, 'z_dof').text = '0'
        self.log.append('added Abaqus BC-VW-mid XSYMM: zero x-displacement on NodeSet BC-VW-mid (328 nodes, x=1.3628)')

    def set_connector_extend(self, extend='constant', names=CONNECTOR_SPRINGS):
        disc = self.root.find('Discrete')
        done = []
        for m in disc.findall('discrete_material'):
            if m.get('name') in names:
                m.find('force').find('extend').text = extend
                done.append(m.get('name'))
        self.log.append(f'connector springs extend={extend}: {done}')

    def add_peb_bottom_to_avw_contact(self, facets):
        """Append the Abaqus Surf-PVW-front-Peb-bottom facets missing from SlidingElastic1Primary."""
        for s in self.mesh.findall('Surface'):
            if s.get('name') == 'SlidingElastic1Primary':
                n0 = len(s)
                for i, c in enumerate(facets, n0 + 1):
                    f = ET.SubElement(s, 'quad4' if len(c) == 4 else 'tri3', {'id': str(i)})
                    f.text = ','.join(str(n) for n in c)
        self.log.append(f'added {len(facets)} PeB-bottom facets to SlidingElastic1Primary (Abaqus Int-AVW-PVW-PEB-inner)')

    def set_diagnostic_output(self):
        ctrl = self.root.find('Control')
        ctrl.find('plot_level').text = 'PLOT_MINOR_ITRS'
        pf = self.root.find('Output').find('plotfile')
        for v in list(pf):
            if v.get('type') not in ('displacement', 'relative volume'):
                pf.remove(v)
        self.log.append('diagnostic output: PLOT_MINOR_ITRS, plot vars displacement + relative volume only')


def set_shell_normal_nodal(model, domains, flag=0):
    """shell_normal_nodal on the named ShellDomains (element normals when 0)."""
    done = []
    for dom in model.root.find('MeshDomains'):
        if dom.tag == 'ShellDomain' and dom.get('name') in domains:
            el = dom.find('shell_normal_nodal')
            if el is None:
                el = ET.SubElement(dom, 'shell_normal_nodal')
            el.text = str(int(flag))
            done.append(dom.get('name'))
    model.log.append(f'shell_normal_nodal={int(flag)} on {done}')
    return done


FANS = ('AVW-Para-L_fan', 'AVW-Para-R_fan', 'CL-L_fan', 'CL-R_fan', 'P-arcus-L_fan', 'P-arcus-R_fan', 'PM_fan',
        'PM_PeB_Left_fan', 'PM_PeB_Right_fan', 'PM_avw_bottom_left_fan', 'PM_avw_bottom_right_fan',
        'PeB-constrin_fan', 'USL-L_fan', 'USL-R_fan')


def remove_stab_springs(model, only_nodeset=None):
    """Delete the zero-length 'stab_*' ground springs (their DiscreteSet, material and binding).
    only_nodeset: restrict to springs whose live node is in this NodeSet (e.g. 'BC-VW-mid')."""
    keep_nodes = None
    if only_nodeset:
        ns = model._nodeset(only_nodeset)
        keep_nodes = {int(v) for v in ns.text.replace('\n', ',').split(',') if v.strip()}
    anch = {int(v) for v in model._nodeset('StabilizationAnchorsSet').text.replace('\n', ',').split(',') if v.strip()}
    names = set()
    for ds in list(model.mesh.findall('DiscreteSet')):
        nm = ds.get('name') or ''
        if not nm.startswith('stab_'):
            continue
        a, b = (int(v) for v in ds.find('delem').text.split(','))
        live = a if b in anch else b
        if keep_nodes is not None and live not in keep_nodes:
            continue
        names.add(nm)
        model.mesh.remove(ds)
    disc = model.root.find('Discrete')
    for el in list(disc):
        if el.get('name') in names or el.get('discrete_set') in names:
            disc.remove(el)
    renumber_discrete_materials(model)
    model.log.append(f'removed {len(names)} stab_* ground springs' + (f' on {only_nodeset}' if only_nodeset else '')
                     + ' (discrete_material ids renumbered 1..N: FEBio reads dmat as a list position)')
    return names


def renumber_discrete_materials(model):
    """FEBio resolves <discrete dmat="k"> as the k-th discrete_material, so ids must stay 1..N in order."""
    disc = model.root.find('Discrete')
    remap = {}
    for i, m in enumerate(disc.findall('discrete_material'), 1):
        remap[m.get('id')] = str(i)
        m.set('id', str(i))
    for b in disc.findall('discrete'):
        b.set('dmat', remap[b.get('dmat')])


def isolate_element_stiffer(model, domain, eid, factor=10.0):
    """DIAGNOSTIC ONLY: move one element into its own domain whose material is the parent's with every
    c_i (and k) multiplied by factor. Tests whether that element drives a failure."""
    eb = model.elem_blocks()
    blk = eb[domain]
    el = next(e for e in blk if int(e.get('id')) == eid)
    blk.remove(el)
    new = ET.Element('Elements', {'type': blk.get('type'), 'name': f'{domain}_e{eid}'})
    new.append(el)
    kids = list(model.mesh)
    model.mesh.insert(kids.index(blk) + 1, new)
    md = model.root.find('MeshDomains')
    parent = next(d for d in md if d.get('name') == domain)
    matname = parent.get('mat')
    mats = model.root.find('Material')
    src = next(m for m in mats if m.get('name') == matname)
    import copy as _copy
    m2 = _copy.deepcopy(src)
    m2.set('name', f'{matname}_x{factor:g}')
    m2.set('id', str(max(int(m.get('id')) for m in mats) + 1))
    for p in m2:
        if p.tag in ('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'k'):
            p.text = '%g' % (float(p.text) * factor)
    mats.append(m2)
    d2 = _copy.deepcopy(parent)
    d2.set('name', f'{domain}_e{eid}')
    d2.set('mat', m2.get('name'))
    md.insert(list(md).index(parent) + 1, d2)
    model.log.append(f'DIAGNOSTIC: element {eid} moved from {domain} to its own domain with {matname} x{factor:g}')
