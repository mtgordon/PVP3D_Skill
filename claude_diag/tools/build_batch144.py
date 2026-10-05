"""Batch 144 (2026-10-02, the user: a new springs-line base with the Yeoh walls, to fit the paper's P1 / P2 from).
From L91_newline_vwyeoh_rhoi0_seamspr001 (faithful Yeoh walls, the canal's sides held by 116 zero-length 0.01 N/mm seam
springs, full wall width and thickness, no tube), add:
  * SlidingElastic1 seg_up 0 -> 2 (solver-side; the answer unchanged on the tube line, healthy and P1);
  * the PM as a structure with its anchors tied to it (build_batch113 step2(INNER + OUTER_BOTTOM)); NOT IN SOURCE;
  * "all interactions on": walls_LA (every outer wall face not in PVW_LA_primary, against the LA shells, as PVW_LA) and
    walls_PM (every outer wall face against the PM shell, PM facets wound toward the walls, offset 1.0); NOT IN SOURCE;
    the new SurfacePairs inserted just before the PVW_LA pair (after the surfaces they name);
  * all 27 PM_conn-family springs on the PM_Plane's inner arc with the mean-scaled law (build_batch129, 129b); NOT IN SOURCE.
  L144_seamspr_pm_allcontact   built, test-read, then run
usage: py -3.10 build_batch144.py"""
import copy
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS, emit  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM  # noqa: E402
from build_batch132 import seg_up2  # noqa: E402
from build_batch131 import facets_xml  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402
import build_batch133 as b133  # noqa: E402
import build_batch136 as b136  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
TOOLS = os.path.dirname(os.path.abspath(__file__))
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
NAME = 'L144_seamspr_pm_allcontact'
PRE, MID = NAME + '_tmppre', NAME + '_tmpmid'
WALLS = (tube.AVW, tube.CX, tube.PVW)


def wall_contacts(model):
    """walls_LA + walls_PM on a model without the tube's wrap; returns the wall nodes that start inside the LA's thickness"""
    tmp = os.path.join(SCRATCH, NAME + '_pre_contact.feb')
    model.write(tmp)
    f, t = Feb(tmp), open(tmp, encoding='utf-8').read()

    def surf(name):
        i = t.index(f'<Surface name="{name}">')
        blk = t[i:t.index('</Surface>', i)]
        return [tuple(int(v) for v in m.split(',')) for m in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', blk)]
    lumen = {tuple(sorted(fc)) for s in ('SlidingElastic1Primary', 'SlidingElastic1Secondary') for fc in surf(s)}
    pvwla = {tuple(sorted(fc)) for fc in surf('PVW_LA_primary')}
    wn = {n for dm in WALLS for n in f.domain_nodes(dm)}
    outer = [fc for fc in tube.boundary_faces(f, WALLS + (tube.PEB,))
             if all(n in wn for n in fc) and tuple(sorted(fc)) not in lumen]
    la_faces = [fc for fc in outer if tuple(sorted(fc)) not in pvwla]
    inside = [n for n in {n for fc in outer for n in fc}
              if (lambda r: r[0] < b133.FAR and r[1] < 2.0)(b133.signed_la(np.array(f.nodes[n])))]
    dom_of = {}
    for dm in WALLS:
        for n in f.domain_nodes(dm):
            dom_of.setdefault(n, dm)
    pm_faces = [tuple(reversed(c)) for c in f.elem_blocks['PM_Plane'][1].values()]
    mesh = model.mesh
    anchor = next(x for x in mesh if x.tag == 'SurfacePair' and x.get('name') == 'PVW_LA')
    for xml in (facets_xml('walls_LA_primary', la_faces), facets_xml('walls_PM_primary', outer),
                facets_xml('PM_contact', pm_faces)):
        mesh.insert(list(mesh).index(anchor), ET.fromstring(xml))
    for nm, sec in (('walls_LA', 'PVW_LA_secondary'), ('walls_PM', 'PM_contact')):
        sp = ET.Element('SurfacePair', name=nm)
        ET.SubElement(sp, 'primary').text = nm + '_primary'
        ET.SubElement(sp, 'secondary').text = sec
        mesh.insert(list(mesh).index(anchor), sp)
    tmpl = next(x for x in model.root.find('Contact') if x.get('name') == 'PVW_LA')
    assert tmpl.findtext('offset').strip() == '2.0'
    for nm, off in (('walls_LA', '2.0'), ('walls_PM', '1.0')):
        con = copy.deepcopy(tmpl)
        con.set('name', nm)
        con.set('surface_pair', nm)
        con.find('offset').text = off
        model.root.find('Contact').append(con)
    model.log.append(f'NOT IN SOURCE: contact walls_LA, {len(la_faces)} outer wall faces not in PVW_LA_primary '
                     f'({dict(Counter(dom_of[fc[0]] for fc in la_faces))}) against the LA shells (PVW_LA_secondary), as '
                     f'PVW_LA (offset 2.0); contact walls_PM, all {len(outer)} outer wall faces against the PM shell '
                     f'(PM_contact: {len(pm_faces)} PM_Plane facets wound toward the walls), as PVW_LA but offset 1.0; '
                     f'wall nodes starting inside the LA\'s thickness: {len(inside)}')
    return inside


if __name__ == '__main__':
    os.makedirs(SCRATCH, exist_ok=True)
    for nm in (NAME, PRE, MID):
        assert not os.path.exists(os.path.join(RUNS, nm)), f'{nm} exists; not overwriting'
    orig = os.path.join(RUNS, BASE, BASE + '.feb')
    m = Model113(orig)
    seg_up2(m)
    step2(INNER + OUTER_BOTTOM)(m)
    inside = wall_contacts(m)
    print(m.log[-1])
    emit(PRE, m)
    b136.run_main(os.path.join(TOOLS, 'build_batch129.py'), {'BASE': PRE, 'NAME': MID})
    b136.run_main(os.path.join(TOOLS, 'build_batch129b.py'), {'SRC': MID, 'OLD': PRE, 'NAME': NAME})
    out = os.path.join(RUNS, NAME)
    notes = [open(os.path.join(RUNS, nm, nm + '.feb.changes.txt'), encoding='utf-8').read() for nm in (PRE, MID)]
    notes.append(open(os.path.join(out, NAME + '.feb.changes.txt'), encoding='utf-8').read())
    with open(os.path.join(out, NAME + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
        fo.write(f'base: {orig} (tools/build_batch144.py: seg_up 2, the PM structure and the wall contacts as {PRE}, then '
                 f'the 27 PM_conn-family springs on the inner arc as {MID}, then their law scaled)\n' + '\n'.join(notes))
    for nm in (PRE, MID):
        shutil.move(os.path.join(RUNS, nm), os.path.join(SCRATCH, nm))
    print('wall nodes inside the LA at rest', len(inside))
    print('wrote', os.path.join(out, NAME + '.feb'))
