"""Batch 134 (2026-10-01, the user: "Does it have a contact between all of the vaginal wall (avw, pvw, cervix, and the
curved sides) and the Levator? If not, add that contact as well"). L133_tube_r30s_pm_narrow (the springs base with the
wrap narrowed clear of the levator, its halves in Vagina_AVW / Vagina_PVW) + NOT IN SOURCE contact walls_LA: the
vaginal walls' outer faces (AVW, cervix, PVW and both wrap halves; boundary faces of the walls, not the lumen faces of
SlidingElastic1, not the faces shared with the perineal body) that are not already in the source's PVW_LA_primary,
against the levator shells (PVW_LA_secondary, all 1170 LA facets). Settings as PVW_LA (sliding-elastic, penalty 0.5 auto,
two_pass, offset 2.0 = half the LA's 4 mm, seg_up 2). With the narrowed wrap no wall face starts inside the LA, so none
is left out (unlike L131_tube_pm_wallcontact, which dropped 274 wrap faces).
  L134_tube_r30s_pm_narrow_wallsLA   built only
usage: py -3.10 build_batch134.py"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402
from build_batch131 import facets_xml  # noqa: E402
import build_batch133 as b133  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE, NAME = 'L133_tube_r30s_pm_narrow', 'L134_tube_r30s_pm_narrow_wallsLA'
WALLS = (tube.AVW, tube.CX, tube.PVW, 'canal_wrap_hex_avw', 'canal_wrap_hex_pvw')

if __name__ == '__main__':
    src = os.path.join(RUNS, BASE, BASE + '.feb')
    d = os.path.join(RUNS, NAME)
    assert not os.path.exists(d), f'{NAME} exists; not overwriting'
    f, t = Feb(src), open(src, encoding='utf-8').read()

    def surf(name):
        i = t.index(f'<Surface name="{name}">')
        blk = t[i:t.index('</Surface>', i)]
        return [tuple(int(v) for v in m.split(',')) for m in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', blk)]
    lumen = {tuple(sorted(fc)) for s in ('SlidingElastic1Primary', 'SlidingElastic1Secondary') for fc in surf(s)}
    pvwla = {tuple(sorted(fc)) for fc in surf('PVW_LA_primary')}
    wn = {n for dm in WALLS for n in f.domain_nodes(dm)}
    faces = [fc for fc in tube.boundary_faces(f, WALLS + (tube.PEB,))
             if all(n in wn for n in fc) and tuple(sorted(fc)) not in lumen and tuple(sorted(fc)) not in pvwla]
    # check: no face node starts inside the LA's thickness
    inside = [n for n in {n for fc in faces for n in fc} if (lambda r: r[0] < b133.FAR and r[1] < 2.0)(b133.signed_la(np.array(f.nodes[n])))]
    assert not inside, f'{len(inside)} wall nodes start inside the LA'
    dom_of = {}
    for dm in WALLS:
        for n in f.domain_nodes(dm):
            dom_of.setdefault(n, dm)
    i = t.index('<contact name="PVW_LA"')
    tmpl = t[i:t.index('</contact>', i) + len('</contact>')]
    con = tmpl.replace('name="PVW_LA" surface_pair="PVW_LA"', 'name="walls_LA" surface_pair="walls_LA"')
    assert con != tmpl
    pair = ('\t\t<SurfacePair name="walls_LA">\n\t\t\t<primary>walls_LA_primary</primary>\n\t\t\t<secondary>PVW_LA_secondary'
            '</secondary>\n\t\t</SurfacePair>\n')
    j = t.index('<SurfacePair name="PVW_LA">')
    j = t.rindex('\n', 0, j) + 1
    t2 = t[:j] + facets_xml('walls_LA_primary', faces) + pair + t[j:]
    k = t2.index('</Contact>')
    t2 = t2[:k] + '\t' + con + '\n\t' + t2[k:]
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
    cnt = Counter(dom_of[fc[0]] for fc in faces)
    msg = (f'{NAME}: {BASE} + NOT IN SOURCE: contact walls_LA, the walls\' outer faces not in PVW_LA_primary ({len(faces)} '
           f'faces: {dict(cnt)}; none starts inside the LA) against the LA shells (PVW_LA_secondary), as PVW_LA otherwise '
           f'(sliding-elastic, penalty 0.5 auto, two_pass, offset 2.0, seg_up 2)\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
