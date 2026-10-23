"""Batch 297 (2026-10-16 evening; the user's loft-contacts plan, step 2: AVW-Para <-> every vaginal wall part (AVW, PVW, cervix,
PeB, curves, apex_contour); memory project-loft-contacts-plan; its own variant, not on top of batch 296).
NOT IN SOURCE: a sliding-elastic contact between both AVW-Para loft sheets (primary) and the walls' outer surface (secondary).
In L294_junction_chains_MAX free sheet nodes go up to 1.2 mm into the curves (canal_wrap_hex_avw / _pvw_fine) by t 0.8
(loft_contacts/para_wall_geom.py); at rest the nearest free sheet nodes are 0.15 mm outside the walls.
  primary: sheet facets with no node on the wall surface (the attached row left out: it sits on the wall), each wound to face
           its nearest wall face (shell contact engages only when the two contact surfaces face each other; gotcha 32)
  secondary: wall boundary faces (outward hex8 faces, all wall domains together) within DMAX mm of a sheet node at rest whose
           nearest sheet node is OUTSIDE them (so the lumen faces, behind the wall, are not in the surface)
Settings as walls_LA (penalty 0.5 auto, seg_up 2), but offset 0, two_pass 0, search_radius SR (default 0.5: below the wall
thickness, so it cannot grab the wall's far side; gotcha 30).
  --all    (L297b) secondary = EVERY wall boundary face within DMAX mm of a free sheet node at rest, except the lumen faces
           (SlidingElastic1Primary / Secondary). L297 had only faces with that node outside them, and so missed the curve-end
           faces near the PeB that face away at rest and turn toward the sheet as the curve end rolls (the nodes that crossed).
  --both   (L297c) also a second contact AVW_Para_walls_b with the primary facets reversed (the sheet's other face): shell
           contact is one-sided (gotcha 32) and some crossing nodes had no primary facet facing the face they crossed.
usage: py -3.10 build_batch297.py SRC_RUN NEW_RUN [--sr 0.5] [--pen 0.5] [--dmax 10] [--all] [--both]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

WALLS = ('_PickedSet347', '_PickedSet66', '_PickedSet64', '_PickedSet346', 'canal_wrap_hex_avw', 'canal_wrap_hex_avw_cervix',
         'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine', 'apex_contour')
SHEETS = ('AVW-Para-L_sheet', 'AVW-Para-R_sheet')
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def main(src, new, sr, pen, dmax, allf=False, both=False):
    t, X, E = load(src)
    cnt, own = {}, {}
    for nm, c in E.values():
        if nm in WALLS:
            for h in HF:
                q = tuple(c[i] for i in h); k = frozenset(q)
                cnt[k] = cnt.get(k, 0) + 1; own[k] = q
    faces = [own[k] for k, v in cnt.items() if v == 1]
    wn = {n for q in faces for n in q}
    FC = np.array([np.mean([X[v] for v in q], 0) for q in faces])
    FN = np.array([np.cross(X[q[2]] - X[q[0]], X[q[3]] - X[q[1]]) for q in faces]); FN /= np.linalg.norm(FN, axis=1)[:, None]
    sheet = [c for nm, c in E.values() if nm in SHEETS]
    sn = sorted({n for c in sheet for n in c if n not in wn})
    SP = np.array([X[n] for n in sn])
    lumen = set()
    if allf:
        import re
        for nm in ('SlidingElastic1Primary', 'SlidingElastic1Secondary'):
            a = t.index(f'<Surface name="{nm}"'); b = t.index('</Surface>', a)
            lumen |= {frozenset(int(v) for v in r.split(',')) for r in re.findall(r'>([0-9,]+)</', t[a:b])}
    sec = []
    for q, cen, nv in zip(faces, FC, FN):
        e = np.linalg.norm(SP - cen, axis=1); j = int(np.argmin(e))
        if e[j] < dmax and (np.dot(SP[j] - cen, nv) > 0 if not allf else frozenset(q) not in lumen):
            sec.append(q)
    SC = np.array([np.mean([X[v] for v in q], 0) for q in sec])
    SN = np.array([np.cross(X[q[2]] - X[q[0]], X[q[3]] - X[q[1]]) for q in sec]); SN /= np.linalg.norm(SN, axis=1)[:, None]
    prim, flipped, skipped = [], 0, 0
    for c in sheet:
        if any(v in wn for v in c):
            skipped += 1; continue
        cen = np.mean([X[v] for v in c], 0)
        j = int(np.argmin(np.linalg.norm(SC - cen, axis=1)))
        n = np.cross(X[c[1]] - X[c[0]], X[c[2]] - X[c[0]])
        if np.dot(n, SC[j] - cen) < 0:
            c = c[::-1]; flipped += 1
        prim.append(c)
    surf = ('\t\t<Surface name="AVW_Para_walls_primary">\n' +
            ''.join(f'\t\t\t<tri3 id="{k}">{",".join(map(str, c))}</tri3>\n' for k, c in enumerate(prim, 1)) +
            '\t\t</Surface>\n\t\t<Surface name="AVW_Para_walls_secondary">\n' +
            ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q))}</quad4>\n' for k, q in enumerate(sec, 1)) +
            '\t\t</Surface>\n')
    if both:
        surf += ('\t\t<Surface name="AVW_Para_walls_primary_b">\n' +
                 ''.join(f'\t\t\t<tri3 id="{k}">{",".join(map(str, c[::-1]))}</tri3>\n' for k, c in enumerate(prim, 1)) +
                 '\t\t</Surface>\n')
    pair = ('\t\t<SurfacePair name="AVW_Para_walls">\n\t\t\t<primary>AVW_Para_walls_primary</primary>\n'
            '\t\t\t<secondary>AVW_Para_walls_secondary</secondary>\n\t\t</SurfacePair>\n')
    if both:
        pair += ('\t\t<SurfacePair name="AVW_Para_walls_b">\n\t\t\t<primary>AVW_Para_walls_primary_b</primary>\n'
                 '\t\t\t<secondary>AVW_Para_walls_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + pair + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', 'name="AVW_Para_walls" surface_pair="AVW_Para_walls"')
         .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>')
         .replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
         .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>')
         .replace('<offset>2.0</offset>', '<offset>0</offset>'))
    for s in ('AVW_Para_walls', f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>',
              '<offset>0</offset>'):
        assert s in c, s
    if both:
        c = c + '\t\t' + c.replace('name="AVW_Para_walls" surface_pair="AVW_Para_walls"',
                                   'name="AVW_Para_walls_b" surface_pair="AVW_Para_walls_b"')
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE: contact AVW_Para_walls (sliding-elastic PENALTY {pen:g} auto, two_pass 0, search_radius {sr:g}, '
           f'offset 0, seg_up 2): primary = {len(prim)} facets of {" + ".join(SHEETS)} ({skipped} facets on the attached row '
           f'left out; {flipped} rewound to face the walls), secondary = {len(sec)} outer wall faces (of {len(faces)} wall '
           f'boundary faces; within {dmax:g} mm of a free sheet node at rest, ' +
           ('all but the lumen faces (SlidingElastic1Primary / Secondary)' if allf else 'that node outside them') + ')' +
           ("; plus contact AVW_Para_walls_b: the same with the primary facets reversed (the sheets' other face)" if both else ''))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch297.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--sr', 0.5), g('--pen', 0.5), g('--dmax', 10.0), '--all' in a, '--both' in a)
