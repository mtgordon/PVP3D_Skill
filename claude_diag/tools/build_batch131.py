"""Batch 131 (2026-10-01, the user: "add contact between all parts of the vaginal wall and the PM_mid and the levator").
On the base L128_tube_r30s_pm (the tube, faithful walls, seg_up 2, the PM as a structure), NOT IN SOURCE (the source's only
wall contact besides the canal is PVW_LA):
  * walls_LA: the vaginal walls' outer faces (AVW, cervix, PVW, the tube's wrap; boundary faces of the walls, not the lumen
    faces of SlidingElastic1, not the faces shared with the perineal body) that are not already in PVW_LA_primary, against
    the levator shells (the existing surface PVW_LA_secondary, all 1170 LA facets). Settings as PVW_LA (sliding-elastic,
    penalty 0.5 auto, two_pass, offset 2.0 = half the LA's 4 mm, seg_up 2). Wrap faces that start inside the LA's
    thickness (a node within 2.5 mm of an LA facet centre, or behind an LA facet within 8 mm) are left out: the wrap
    (NOT IN SOURCE itself) overlaps the LA at rest there, and a contact would push it out at t = 0.
  * walls_PM: all the walls' outer faces against the PM shell (PM_Plane, 2 mm; offset 1.0), same settings. The walls lie
    on the PM's back side, so the PM facets are written with reversed winding (their normal toward the walls).
  L131_tube_pm_wallcontact   L128_tube_r30s_pm + walls_LA + walls_PM
usage: py -3.10 build_batch131.py"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE, NAME = 'L128_tube_r30s_pm', 'L131_tube_pm_wallcontact'
WALLS = (tube.AVW, tube.CX, tube.PVW, 'canal_wrap_hex')
LA = ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri')


def facets_xml(name, faces):
    out = [f'\t\t<Surface name="{name}">']
    for i, fc in enumerate(faces, 1):
        tag = 'quad4' if len(fc) == 4 else 'tri3'
        out.append(f'\t\t\t<{tag} id="{i}">{",".join(map(str, fc))}</{tag}>')
    out.append('\t\t</Surface>')
    return '\n'.join(out) + '\n'


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
    wrapn = set(f.domain_nodes('canal_wrap_hex'))
    faces = [fc for fc in tube.boundary_faces(f, WALLS + (tube.PEB,))
             if all(n in wn for n in fc) and tuple(sorted(fc)) not in lumen]
    # LA facets: centres and normals
    la = [c for dm in LA for c in f.elem_blocks[dm][1].values()]
    C = np.array([np.mean([f.nodes[n] for n in c], 0) for c in la])
    P = np.array([[f.nodes[n] for n in c[:3]] for c in la])
    N = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    N /= np.linalg.norm(N, axis=1)[:, None]

    def inside_la(n):
        x = np.array(f.nodes[n])
        k = int(np.argmin(np.linalg.norm(C - x, axis=1)))
        dc = np.linalg.norm(C[k] - x)
        return dc < 2.5 or (dc < 8 and (x - C[k]) @ N[k] < 0)
    bad = {n for n in {n for fc in faces for n in fc} & wrapn if inside_la(n)}
    la_faces = [fc for fc in faces if tuple(sorted(fc)) not in pvwla and not any(n in bad for n in fc)]
    dropped = [fc for fc in faces if tuple(sorted(fc)) not in pvwla and any(n in bad for n in fc)]
    dom_of = {}
    for dm in WALLS:
        for n in f.domain_nodes(dm):
            dom_of.setdefault(n, dm)
    pm_faces = [tuple(reversed(c)) for c in f.elem_blocks['PM_Plane'][1].values()]
    contact = open(src, encoding='utf-8').read()
    i = t.index('<contact name="PVW_LA"')
    tmpl = t[i:t.index('</contact>', i) + len('</contact>')]

    def make(name, offset):
        c = tmpl.replace('name="PVW_LA" surface_pair="PVW_LA"', f'name="{name}" surface_pair="{name}"')
        c = c.replace('<offset>2.0</offset>', f'<offset>{offset}</offset>')
        assert c != tmpl and f'name="{name}"' in c
        return c
    surfs = facets_xml('walls_LA_primary', la_faces) + facets_xml('walls_PM_primary', faces) + facets_xml('PM_contact', pm_faces)
    pairs = ('\t\t<SurfacePair name="walls_LA">\n\t\t\t<primary>walls_LA_primary</primary>\n\t\t\t<secondary>PVW_LA_secondary'
             '</secondary>\n\t\t</SurfacePair>\n\t\t<SurfacePair name="walls_PM">\n\t\t\t<primary>walls_PM_primary</primary>'
             '\n\t\t\t<secondary>PM_contact</secondary>\n\t\t</SurfacePair>\n')
    j = t.index('<SurfacePair name="PVW_LA">')
    j = t.rindex('\n', 0, j) + 1
    t2 = t[:j] + surfs + pairs + t[j:]
    k = t2.index('</Contact>')
    t2 = t2[:k] + make('walls_LA', '2.0') + '\n\t\t' + make('walls_PM', '1.0') + '\n\t' + t2[k:]
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
    cnt = Counter(dom_of[fc[0]] for fc in la_faces)
    msg = (f'{NAME}: {BASE} + NOT IN SOURCE: contact walls_LA (the walls\' outer faces not in PVW_LA_primary: {len(la_faces)} '
           f'faces, {dict(cnt)}; {len(dropped)} wrap faces that start inside the LA\'s thickness left out) against the LA '
           f'shells (PVW_LA_secondary), offset 2.0; contact walls_PM (all {len(faces)} wall outer faces) against the PM shell '
           f'(PM_contact: {len(pm_faces)} PM_Plane facets, winding reversed to face the walls), offset 1.0; both as PVW_LA '
           f'otherwise (sliding-elastic, penalty 0.5 auto, two_pass, seg_up 2)\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
