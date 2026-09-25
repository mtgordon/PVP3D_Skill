"""Build single-variable test variants of the PVP3D FEBio model.

Every operation edits an ElementTree of the .feb and keeps FEBio's ordering rules:
new nodes go in a new <Nodes> block right after the block holding the current max
node id, and anything referencing new nodes (DiscreteSet) is appended at the end
of <Mesh>.
"""
import copy
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from inpmap import part_nodes, nearest_feb  # noqa: E402
from febmodel import Feb  # noqa: E402
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import JOBS_DIR

JOBS = JOBS_DIR

# Abaqus LA-new -> VW-PeB axial connectors (Sphincter-L-1..4, Sphincter-R-1..4)
SPHINCTER_CONN = [(3, 132), (5, 134), (7, 136), (19, 129), (24, 155), (26, 153), (28, 151), (39, 149)]
# *Connector Elasticity, nonlinear, LA-Y-sphincter-side-*: (force, displacement)
SPHINCTER_CURVE = [(0, 0), (0.7, 4.7), (1.68, 9.515), (3.36, 13.985), (5.04, 18.8), (7.56, 23.5),
                   (10.5, 28.2), (14.56, 33.015), (19.6, 37.715), (26.46, 42.415), (35.56, 47.115)]


class Model:
    def __init__(self, path):
        self.src = path
        self.tree = ET.parse(path)
        self.root = self.tree.getroot()
        self.mesh = self.root.find('Mesh')
        self.log = []

    # ---------- helpers ----------
    def nodes(self):
        out = {}
        for blk in self.mesh.findall('Nodes'):
            for n in blk:
                out[int(n.get('id'))] = np.array([float(v) for v in n.text.split(',')])
        return out

    def max_node_id(self):
        return max(int(n.get('id')) for blk in self.mesh.findall('Nodes') for n in blk)

    def elem_blocks(self):
        return {e.get('name'): e for e in self.mesh.findall('Elements')}

    def add_nodes(self, name, coords):
        """Append new nodes (list of xyz) after the max-id Nodes block; return new ids."""
        blocks = self.mesh.findall('Nodes')
        mx_blk = max(blocks, key=lambda b: max(int(n.get('id')) for n in b))
        start = self.max_node_id() + 1
        new = ET.Element('Nodes', {'name': name})
        new.text = '\n\t\t\t'
        ids = []
        for i, xyz in enumerate(coords):
            e = ET.SubElement(new, 'node', {'id': str(start + i)})
            e.text = ','.join('%.9g' % v for v in xyz)
            e.tail = '\n\t\t\t'
            ids.append(start + i)
        new.tail = '\n\t\t'
        idx = list(self.mesh).index(mx_blk)
        self.mesh.insert(idx + 1, new)
        return ids

    def add_discrete_set(self, name, pairs):
        ds = ET.SubElement(self.mesh, 'DiscreteSet', {'name': name})
        ds.text = '\n\t\t\t'
        for a, b in pairs:
            d = ET.SubElement(ds, 'delem')
            d.text = f'{a},{b}'
            d.tail = '\n\t\t\t'
        ds.tail = '\n\t'

    def next_dmat_id(self):
        disc = self.root.find('Discrete')
        return max(int(m.get('id')) for m in disc.findall('discrete_material')) + 1

    def add_nonlinear_spring(self, name, set_name, pts, measure='elongation', extend='extrapolate'):
        disc = self.root.find('Discrete')
        mid = self.next_dmat_id()
        mats = disc.findall('discrete_material')
        m = ET.Element('discrete_material', {'id': str(mid), 'name': name, 'type': 'nonlinear spring'})
        ET.SubElement(m, 'scale').text = '1'
        ET.SubElement(m, 'measure').text = measure
        fo = ET.SubElement(m, 'force', {'type': 'point'})
        ET.SubElement(fo, 'interpolate').text = 'linear'
        ET.SubElement(fo, 'extend').text = extend
        p = ET.SubElement(fo, 'points')
        for x, y in pts:
            ET.SubElement(p, 'pt').text = f'{x},{y}'
        disc.insert(list(disc).index(mats[-1]) + 1, m)
        b = ET.SubElement(disc, 'discrete', {'dmat': str(mid), 'discrete_set': set_name})
        b.tail = '\n\t'
        return mid

    def add_linear_spring(self, name, set_name, E):
        disc = self.root.find('Discrete')
        mid = self.next_dmat_id()
        mats = disc.findall('discrete_material')
        m = ET.Element('discrete_material', {'id': str(mid), 'name': name, 'type': 'linear spring'})
        ET.SubElement(m, 'E').text = str(E)
        disc.insert(list(disc).index(mats[-1]) + 1, m)
        ET.SubElement(disc, 'discrete', {'dmat': str(mid), 'discrete_set': set_name})
        return mid

    # ---------- operations ----------
    def set_control(self, **kw):
        ctrl = self.root.find('Control')
        for k, v in kw.items():
            el = ctrl.find('.//' + k)
            assert el is not None, k
            el.text = str(v)
        self.log.append(f'control {kw}')

    def set_contact(self, name_substr, **kw):
        n = 0
        for c in self.root.find('Contact'):
            if name_substr in c.get('name'):
                for k, v in kw.items():
                    el = c.find(k)
                    if v is None:
                        if el is not None:
                            c.remove(el)
                        continue
                    if el is None:
                        el = ET.SubElement(c, k)
                    if isinstance(v, tuple):     # (value, attrib dict)
                        el.text = str(v[0])
                        el.attrib.clear(); el.attrib.update(v[1])
                    else:
                        el.text = str(v)
                        if k == 'penalty':
                            el.attrib.pop('lc', None)
                n += 1
        self.log.append(f'contact *{name_substr}* x{n}: {kw}')
        return n

    def remove_contact(self, name_substr):
        con = self.root.find('Contact')
        for c in list(con):
            if name_substr in c.get('name'):
                con.remove(c)
                self.log.append(f'removed contact {c.get("name")}')

    def unmerge_la(self, other_domains=('Sphincter-L_fan', 'Sphincter-R_fan')):
        """Give the LA its own copies of nodes it shares with other_domains."""
        eb = self.elem_blocks()
        la_names = [n for n in eb if n.startswith('LA_')]
        la_nodes = set()
        for n in la_names:
            for e in eb[n]:
                la_nodes.update(int(v) for v in e.text.split(','))
        oth = set()
        for n in other_domains:
            for e in eb[n]:
                oth.update(int(v) for v in e.text.split(','))
        shared = sorted(la_nodes & oth)
        X = self.nodes()
        new_ids = self.add_nodes('LA_unmerged_nodes', [X[i] for i in shared])
        remap = dict(zip(shared, new_ids))
        for n in la_names:
            for e in eb[n]:
                e.text = ','.join(str(remap.get(int(v), int(v))) for v in e.text.split(','))
        # LA-only node sets / surfaces follow the LA copy
        for ns in self.mesh.findall('NodeSet'):
            if ns.get('name').startswith('BC-LA'):
                ids = [int(v) for v in ns.text.replace('\n', ',').split(',') if v.strip()]
                if any(i in remap for i in ids):
                    ns.text = ','.join(str(remap.get(i, i)) for i in ids)
        for s in self.mesh.findall('Surface'):
            if '_LA_ICM_tie_primary' in s.get('name'):
                for f in s:
                    f.text = ','.join(str(remap.get(int(v), int(v))) for v in f.text.split(','))
        self.log.append(f'unmerged {len(shared)} LA nodes from {other_domains}: {shared} -> {new_ids}')
        self.remap_la = remap
        return remap

    def add_sphincter_connectors(self, curve=SPHINCTER_CURVE, scale=1.0, extend='extrapolate'):
        """Abaqus LA-Y-sphincter-side connectors as FEBio nonlinear springs (LA node -> PeB node)."""
        la = part_nodes('LA-new')
        vw = part_nodes('VW-PeB')
        X = self.nodes()
        remap = getattr(self, 'remap_la', {})
        pairs = []
        for a, b in SPHINCTER_CONN:
            fa, da = nearest_feb_dict(X, la[a])
            fb, db = nearest_feb_dict(X, vw[b])
            assert da < 1e-3 and db < 1e-3, (a, b, da, db)
            pairs.append((remap.get(fa, fa), fb))
        self.add_discrete_set('LA_sphincter_connectors', pairs)
        pts = [(u, F * scale) for F, u in curve]
        self.add_nonlinear_spring('LA_sphincter_conn_mat', 'LA_sphincter_connectors', pts, extend=extend)
        self.log.append(f'added {len(pairs)} Abaqus sphincter connectors {pairs} scale={scale}')
        return pairs

    def remove_truss_ties(self):
        """Strip the phantom-ribbon tied contacts and phantom links from the TRUSS model."""
        con = self.root.find('Contact')
        names = [c.get('name') for c in con if '_LA_ICM_tie' in c.get('name')]
        for c in list(con):
            if c.get('name') in names:
                con.remove(c)
        for el in list(self.mesh):
            nm = el.get('name') or ''
            if (el.tag in ('Surface', 'SurfacePair') and '_LA_ICM_tie' in nm) or \
               (el.tag == 'DiscreteSet' and nm.endswith('_phantom_links')):
                self.mesh.remove(el)
        disc = self.root.find('Discrete')
        for el in list(disc):
            if (el.get('discrete_set') or '').endswith('_phantom_links') or el.get('name') == 'tie_phantom_link_mat':
                disc.remove(el)
        self.log.append(f'removed tied contacts {names}, phantom links + their material')

    def add_abaqus_ties_lc(self, penalty=10, tol=0.01, maxaug=10):
        """Abaqus *Tie LA-new (slave) -> ATLA/Posterior_Arcus truss nodes (master) as exact
        FEBio linear constraints u_LA = u_truss (x, y, z)."""
        from abq_ties import resolve
        res, la = resolve()
        X = self.nodes()
        eb = self.elem_blocks()
        la_ids = sorted({int(v) for n in eb if n.startswith('LA_') for e in eb[n] for v in e.text.split(',')})
        LAX = {i: X[i] for i in la_ids}
        truss = {int(n.get('id')): X[int(n.get('id'))] for blk in self.mesh.findall('Nodes')
                 if (blk.get('name') or '').endswith('_centerline') for n in blk}
        best = {}   # slave feb id -> (dist, master feb id, tie)
        for name, (pairs, mn) in res.items():
            for s, m, d in pairs:
                fs, ds = nearest_feb_dict(LAX, la[s])
                fm, dm = nearest_feb_dict(truss, mn[m])
                assert ds < 1e-3 and dm < 1e-3, (name, s, m, ds, dm)
                if fs not in best or d < best[fs][0]:
                    best[fs] = (d, fm, name)
        cons = ET.Element('Constraints')
        c = ET.SubElement(cons, 'constraint', {'name': 'LA_truss_ties', 'type': 'linear constraint'})
        ET.SubElement(c, 'tol').text = str(tol)
        ET.SubElement(c, 'penalty').text = str(penalty)
        ET.SubElement(c, 'maxaug').text = str(maxaug)
        for fs, (d, fm, name) in sorted(best.items()):
            for dof in 'xyz':
                lc = ET.SubElement(c, 'linear_constraint')
                ET.SubElement(lc, 'node', {'id': str(fs), 'bc': dof}).text = '1'
                ET.SubElement(lc, 'node', {'id': str(fm), 'bc': dof}).text = '-1'
        kids = list(self.root)
        after = self.root.find('Contact')
        self.root.insert(kids.index(after) + 1, cons)
        per = {}
        for fs, (d, fm, name) in best.items():
            per[name] = per.get(name, 0) + 1
        self.log.append(f'added Abaqus LA-truss ties as linear constraints (penalty={penalty}, tol={tol}, '
                        f'maxaug={maxaug}): {per}')
        return best

    def add_all_sphincter_connectors(self, which=('L', 'R', 'P'), scale=1.0):
        """All Abaqus LA-new -> VW-PeB connectors: Sphincter-L-1..8, -R-1..8, -P-1..24."""
        conns = {
            'L': [(3, 132), (5, 134), (7, 136), (19, 129), (2, 131), (4, 133), (6, 135), (20, 137)],
            'R': [(24, 155), (26, 153), (28, 151), (39, 149), (23, 156), (25, 154), (27, 152), (40, 150)],
            'P': [(18, 128), (17, 127), (16, 126), (15, 125), (14, 124), (13, 123), (12, 122), (11, 121),
                  (10, 120), (9, 119), (44, 118), (43, 117), (43, 157), (42, 158), (29, 159), (30, 160),
                  (31, 161), (32, 162), (33, 163), (34, 164), (35, 165), (36, 166), (37, 167), (38, 168)],
        }
        # *Connector Elasticity (force, displacement) for each behavior family
        side = SPHINCTER_CURVE
        post = [(0, 0), (11.6667, 4.7), (28, 9.515), (56, 13.985), (84, 18.8), (126, 23.5), (175, 28.2),
                (242.667, 33.015), (326.667, 37.715), (441, 42.415), (592.667, 47.115)]
        la = part_nodes('LA-new')
        vw = part_nodes('VW-PeB')
        X = self.nodes()
        eb = self.elem_blocks()
        la_ids = sorted({int(v) for n in eb if n.startswith('LA_') for e in eb[n] for v in e.text.split(',')})
        pe_ids = sorted({int(v) for n, blk in eb.items() if blk.get('type') == 'hex8'
                         for e in blk for v in e.text.split(',')})
        LAX = {i: X[i] for i in la_ids}
        PEX = {i: X[i] for i in pe_ids}
        for fam, curve in (('side', side), ('post', post)):
            keys = [k for k in which if (k == 'P') == (fam == 'post')]
            pairs = []
            for k in keys:
                for a, b in conns[k]:
                    fa, da = nearest_feb_dict(LAX, la[a])
                    fb, db = nearest_feb_dict(PEX, vw[b])
                    assert da < 1e-3 and db < 1e-3, (k, a, b, da, db)
                    pairs.append((fa, fb))
            if not pairs:
                continue
            nm = f'LA_sphincter_{fam}_conn'
            self.add_discrete_set(nm, pairs)
            self.add_nonlinear_spring(nm + '_mat', nm, [(u, F * scale) for F, u in curve])
            self.log.append(f'added {len(pairs)} Abaqus {fam} sphincter connectors ({"+".join(keys)}), scale={scale}')

    def _nodeset(self, name):
        for ns in self.mesh.findall('NodeSet'):
            if ns.get('name') == name:
                return ns
        raise KeyError(name)

    def _extend_nodeset(self, name, ids):
        ns = self._nodeset(name)
        cur = [int(v) for v in (ns.text or '').replace('\n', ',').split(',') if v.strip()]
        new = [i for i in ids if i not in cur]
        ns.text = ','.join(str(i) for i in cur + new)
        return new

    def restore_remesh_lost_pins(self):
        """Re-attach BCs lost when the v3 hybrid remesh renumbered fan boundary nodes:
        AVW-Para fan apex (BC-RP-ParaSupport) and PeB-constrin_fan anchor edge
        (BC-RP-PeB-constrain-origins). Matches by exact coordinate."""
        X = self.nodes()
        eb = self.elem_blocks()

        def dom_nodes(name):
            return sorted({int(v) for e in eb[name] for v in e.text.split(',')})
        added = {}
        # 1) ParaSupport apex nodes
        ns = self._nodeset('BC-RP-ParaSupport')
        ids = [int(v) for v in ns.text.replace('\n', ',').split(',') if v.strip()]
        cand = dom_nodes('AVW-Para-L_fan') + dom_nodes('AVW-Para-R_fan')
        hits = [c for c in cand for o in ids if np.linalg.norm(X[c] - X[o]) < 1e-5]
        added['BC-RP-ParaSupport'] = self._extend_nodeset('BC-RP-ParaSupport', hits)
        # 2) PeB-constrain anchor edge: fan nodes on the straight origin line
        ns = self._nodeset('BC-RP-PeB-constrain-origins')
        ids = [int(v) for v in ns.text.replace('\n', ',').split(',') if v.strip()]
        O = np.array([X[o] for o in ids])
        a, b = O[np.argmin(O[:, 0])], O[np.argmax(O[:, 0])]
        u = (b - a) / np.linalg.norm(b - a)
        hits = []
        for c in dom_nodes('PeB-constrin_fan'):
            p = X[c] - a
            s = p @ u
            if -1e-5 <= s <= np.linalg.norm(b - a) + 1e-5 and np.linalg.norm(p - s * u) < 1e-5:
                hits.append(c)
        added['BC-RP-PeB-constrain-origins'] = self._extend_nodeset('BC-RP-PeB-constrain-origins', hits)
        self.log.append(f'restored remesh-lost pins: {added}')
        return added

    def remove_domains(self, names):
        """Delete <Elements> + <ShellDomain>/<SolidDomain> for the named domains (nodes stay)."""
        for el in list(self.mesh):
            if el.tag == 'Elements' and el.get('name') in names:
                self.mesh.remove(el)
        md = self.root.find('MeshDomains')
        for el in list(md):
            if el.get('name') in names:
                md.remove(el)
        self.log.append(f'removed domains {names}')

    # Abaqus *Density per material (mm-N-s-tonne), keyed by FEBio material name
    ABQ_DENSITY = {
        'Vagina_AVW': 1.06e-09, 'Vagina_PVW': 1.06e-09, 'Vagina_Cervix': 1.06e-09,
        'PeB-Vagina500%stiffer': 1.06e-09, 'LA_Yamada50pct_Ogden': 1.06e-09,
        'ATLA-Hyper': 0.00011, 'Posterior_Arcus-Hyper': 0.00011, 'Beam-CL-USL': 7.8e-07,
    }

    def set_dynamic(self, rhoi=0.5, mass_scale=1.0):
        """Implicit dynamic analysis with the Abaqus densities (the source model is *Dynamic, Explicit)."""
        ctrl = self.root.find('Control')
        ctrl.find('analysis').text = 'DYNAMIC'
        ctrl.find('.//rhoi').text = str(rhoi)
        for mat in self.root.find('Material'):
            rho = self.ABQ_DENSITY.get(mat.get('name'))
            assert rho is not None, mat.get('name')
            mat.find('density').text = '%g' % (rho * mass_scale)
        self.log.append(f'DYNAMIC analysis, rhoi={rhoi}, Abaqus densities x{mass_scale}')

    def set_la_material(self, **kw):
        for mat in self.root.find('Material'):
            if mat.get('name') == 'LA_Yamada50pct_Ogden':
                for k, v in kw.items():
                    mat.find(k).text = str(v)
        self.log.append(f'LA material {kw}')

    def set_la_shell_type(self, typ):
        for dom in self.root.find('MeshDomains'):
            if dom.get('name', '').startswith('LA_'):
                dom.set('type', typ)
        self.log.append(f'LA ShellDomain type={typ}')

    def set_la_shell_normal_nodal(self, flag):
        for dom in self.root.find('MeshDomains'):
            if dom.get('name', '').startswith('LA_'):
                el = dom.find('shell_normal_nodal')
                if el is None:
                    el = ET.SubElement(dom, 'shell_normal_nodal')
                el.text = str(int(flag))
        self.log.append(f'LA shell_normal_nodal={int(flag)}')

    def write(self, path):
        ET.indent(self.tree, space='\t') if hasattr(ET, 'indent') else None
        self.tree.write(path, encoding='ISO-8859-1', xml_declaration=True)
        with open(path + '.changes.txt', 'w') as fh:
            fh.write(f'base: {self.src}\n')
            for l in self.log:
                fh.write(l + '\n')


def nearest_feb_dict(X, xyz):
    ids = np.fromiter(X.keys(), dtype=int)
    P = np.array([X[i] for i in ids])
    d = np.linalg.norm(P - xyz, axis=1)
    k = int(np.argmin(d))
    return int(ids[k]), float(d[k])


def emit(name, model=None, src=None):
    import shutil
    out = os.path.join(JOBS, 'claude_diag', 'runs')
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    dst = os.path.join(d, name + '.feb')
    if model is None:
        shutil.copy(src, dst)
    else:
        model.write(dst)
    print('wrote', dst)


def set_la_thickness_equiv(model, t_new, t_old=4.0):
    """Thinner LA shell with membrane stiffness preserved (c1, k scaled by t_old/t_new)."""
    s = t_old / t_new
    for dom in model.root.find('MeshDomains'):
        if dom.get('name', '').startswith('LA_'):
            dom.find('shell_thickness').text = str(t_new)
    for mat in model.root.find('Material'):
        if mat.get('name') == 'LA_Yamada50pct_Ogden':
            for k in ('c1', 'k'):
                mat.find(k).text = '%.6g' % (float(mat.find(k).text) * s)
    model.log.append(f'LA thickness {t_old}->{t_new} with c1,k x{s:.4g} (membrane-equivalent)')
