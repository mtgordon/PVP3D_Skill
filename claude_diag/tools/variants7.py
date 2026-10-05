"""Seventh round (2026-09-27, Test A): the user's fitted lofts put back into the springs line, one family at a time.

Why (NEXT_SESSION_PROMPT_2026-09-26b.md, Test A): the lofts line (L26_lofts_la3_pm) runs faster than the springs line
(L26_springs_la3_pm) with a third fewer Newton iterations; which of the four fitted loft families (AVW-Para, CL, USL, PM)
makes the difference?
New operation:
  connectors_to_lofts  the reverse of variants6.Model6.lofts_to_connectors, copied from a model that has the lofts
                       (L26_lofts_la3_pm; both lines' <Nodes> blocks hash the same, so the loft nodes keep their ids).
                       For each family: its <Elements> block and <ShellDomain> (thickness, shell_normal_nodal 0, material
                       *_fan_fit, which the springs line still defines but does not use) and every contact with a surface
                       lying on the loft (the tied-node-on-facet ties of AVW-Para-L and USL-L) with its <SurfacePair> and
                       <Surface>s, each inserted where it sits in the lofts model; the family's springs (DiscreteSet
                       <fam>_conn[_k], its discrete_material and <discrete>) removed and the discrete materials renumbered
                       (FEBio reads dmat="k" as the k-th discrete_material); the CL/USL BC node sets restored to the lofts
                       model's (lofts_to_connectors added the CL/USL line nodes, the springs' anchors, which are in no
                       element in either model).
"""
import copy
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(__file__))
from variants import JOBS, emit  # noqa: E402,F401
from variants2 import renumber_discrete_materials  # noqa: E402
from variants6 import Model6, LINE_BC  # noqa: E402

# family -> (loft domains, connector families as lofts_to_connectors named them)
FAMILIES = {
    'AVW-Para': (('AVW-Para-L_fan', 'AVW-Para-R_fan'), ('AVW-Para-L', 'AVW-Para-R')),
    'CL': (('CL-L_fan', 'CL-R_fan'), ('CL-L', 'CL-R')),
    'USL': (('USL-L_fan', 'USL-R_fan'), ('USL-L', 'USL-R')),
    'PM': (('PM_fan',), ('PM',)),
}


def _key(el):
    return el.tag, el.get('name')


def _insert_like(src_parent, src_el, dst_parent):
    """Insert a copy of src_el into dst_parent after the nearest preceding sibling (in src_parent) that dst_parent also
    has, matched by (tag, name); at the start if there is none."""
    have = {_key(e): i for i, e in enumerate(dst_parent)}
    assert _key(src_el) not in have, f'{_key(src_el)} already in the target'
    sibs = list(src_parent)
    pos = 0
    for prev in reversed(sibs[:sibs.index(src_el)]):
        if _key(prev) in have and prev.get('name') is not None:
            pos = have[_key(prev)] + 1
            break
    dst_parent.insert(pos, copy.deepcopy(src_el))


def _canon(el):
    return re.sub(r'\s+', '', ET.tostring(el, encoding='unicode'))


class Model7(Model6):
    def connectors_to_lofts(self, lofts_feb, family):
        lofts, conn_fams = FAMILIES[family]
        src = Model6(lofts_feb)
        smesh, sdoms, scont = src.mesh, src.root.find('MeshDomains'), src.root.find('Contact')
        # the loft nodes must be the same nodes in both models
        S, T = src.nodes(), self.nodes()
        for loft in lofts:
            blk = src.elem_blocks()[loft]
            ln = {int(v) for e in blk for v in e.text.split(',')}
            assert all(n in T and (abs(T[n] - S[n]) < 1e-12).all() for n in ln), f'{loft}: loft nodes differ'
            assert loft not in self.elem_blocks(), f'{loft} already in the target'
            # elements, shell domain, material (already defined, and identical)
            _insert_like(smesh, blk, self.mesh)
            dom = next(d for d in sdoms if d.get('name') == loft)
            _insert_like(sdoms, dom, self.root.find('MeshDomains'))
            smat = src._material(dom.get('mat'))
            tmat = self._material(dom.get('mat'))
            assert tmat is not None and _canon(tmat) == _canon(smat), f'{loft}: material {dom.get("mat")} differs'
            # contacts with a surface on the loft, with their surface pair and surfaces
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
                moved.append(f'{c.get("name")} ({c.get("type")}, {names[0]} {len(src.surface_facets(names[0]))} '
                             f'facets | {names[1]} {len(src.surface_facets(names[1]))} facets)')
            self.log.append(f'{loft}: loft back from {os.path.basename(lofts_feb)}: {len(blk)} {blk.get("type")} elements '
                            f'{blk[0].get("id")}-{blk[-1].get("id")} on {len(ln)} nodes, ShellDomain mat={dom.get("mat")} '
                            f'thickness {dom.find("shell_thickness").text}, shell_normal_nodal '
                            f'{dom.find("shell_normal_nodal").text}; contacts with a surface on it: {moved or "none"}')
        # the family's springs out
        disc = self.root.find('Discrete')
        dm = disc.findall('discrete_material')
        assert [int(m.get('id')) for m in dm] == list(range(1, len(dm) + 1)), 'dmat ids must be list positions'
        for fam in conn_fams:
            pat = re.compile(re.escape(fam.rstrip('_')) + r'_conn(_\d+)?$')
            sets = [d for d in self.mesh.findall('DiscreteSet') if pat.match(d.get('name'))]
            assert sets, f'{fam}: no spring sets'
            gone = []
            for ds in sets:
                nm = ds.get('name')
                b = next(b for b in disc.findall('discrete') if b.get('discrete_set') == nm)
                # by id: ids stay as written until the renumbering at the end
                mat = next(m for m in disc.findall('discrete_material') if m.get('id') == b.get('dmat'))
                assert mat.get('name') == nm + '_mat', (nm, mat.get('name'))
                disc.remove(b)
                disc.remove(mat)
                self.mesh.remove(ds)
                gone.append(f'{nm} ({len(ds)} springs)')
            self.log.append(f'{fam}: its Abaqus connectors as springs removed: {", ".join(gone)}')
            if fam in LINE_BC:
                bc = LINE_BC[fam]
                old, new = set(self.nodesets_ids(bc)), set(src.nodesets_ids(bc))
                extra = sorted(old - new)
                assert new <= old, (bc, sorted(new - old))
                used = {int(v) for blk in self.elem_blocks().values() for e in blk for v in e.text.split(',')}
                used |= {int(v) for d in self.mesh.findall('DiscreteSet') for e in d for v in e.text.split(',')}
                assert not (set(extra) & used), (bc, sorted(set(extra) & used))
                cur = self._nodeset(bc)
                i = list(self.mesh).index(cur)
                self.mesh.remove(cur)
                self.mesh.insert(i, copy.deepcopy(src._nodeset(bc)))
                self.log.append(f'{bc}: back to the lofts model\'s {len(new)} nodes (the {len(extra)} CL/USL line nodes '
                                f'{extra}, the springs\' anchors, are now in no element or spring)')
        renumber_discrete_materials(self)
        self.log.append(f'discrete materials renumbered 1..{len(disc.findall("discrete_material"))} (dmat = list position)')
