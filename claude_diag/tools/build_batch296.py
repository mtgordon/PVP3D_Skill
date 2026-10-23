"""Batch 296 (2026-10-16 evening; the user's loft-contacts plan, step 1: "AVW-Para <-> PM arch"; memory project-loft-contacts-plan).
NOT IN SOURCE: a sliding-elastic contact between both AVW-Para loft sheets (AVW-Para-L_sheet + AVW-Para-R_sheet, primary) and
the PM arch shell (PM_Plane, secondary). In L294_junction_chains_MAX the sheets pass through the PM arch (at rest every sheet node
near the PM is 0.6-3.2 mm on its -normal side; by t 0.8 22-41 a side are on the + side; loft_contacts/para_pm_geom.py).
Shell contact (shell_tests/gen_ss2.py, skill gotcha 32): it engages when both contact SURFACES face each other; the surface
winding counts, not the element normals. So: each sheet facet is wound to face its nearest PM facet's +normal (toward the PM),
and the PM facets are written reversed (facing the sheets on the -side). shell_bottom_* 0: contact at the PM's node plane
(bottom 1 would put it 2 mm (the PM thickness) toward the sheets, inside the 0.6 mm rest gap).
Settings as walls_LA (penalty 0.5 auto, seg_up 2), but offset 0, two_pass 0, search_radius SR (default 1.0).
usage: py -3.10 build_batch296.py SRC_RUN NEW_RUN [--sr 1.0] [--pen 0.5] [--tp 0]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

SHEETS = ('AVW-Para-L_sheet', 'AVW-Para-R_sheet')


def nrm(X, c):
    a, b, d = (X[v] for v in c[:3])
    n = np.cross(b - a, d - a)
    return n / np.linalg.norm(n)


def main(src, new, sr, pen, tp):
    t, X, E = load(src)
    pm = [c for nm, c in E.values() if nm == 'PM_Plane']
    PC = np.array([np.mean([X[v] for v in c], 0) for c in pm])
    PN = np.array([nrm(X, c) for c in pm])
    prim, flipped, side = [], 0, []
    for nm in SHEETS:
        for c in [c for n2, c in E.values() if n2 == nm]:
            cen = np.mean([X[v] for v in c], 0)
            j = int(np.argmin(np.linalg.norm(PC - cen, axis=1)))
            side.append(np.dot(cen - PC[j], PN[j]))
            if np.dot(nrm(X, c), PN[j]) < 0:
                c = c[::-1]; flipped += 1
            prim.append(c)
    side = np.array(side)
    surf = ('\t\t<Surface name="AVW_Para_PM_primary">\n' +
            ''.join(f'\t\t\t<tri3 id="{k}">{",".join(map(str, c))}</tri3>\n' for k, c in enumerate(prim, 1)) +
            '\t\t</Surface>\n\t\t<Surface name="AVW_Para_PM_secondary">\n' +
            ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, c[::-1]))}</quad4>\n' for k, c in enumerate(pm, 1)) +
            '\t\t</Surface>\n')
    pair = ('\t\t<SurfacePair name="AVW_Para_PM">\n\t\t\t<primary>AVW_Para_PM_primary</primary>\n'
            '\t\t\t<secondary>AVW_Para_PM_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + pair + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    blk = t[a:b]
    c = (blk.replace('name="walls_LA" surface_pair="walls_LA"', 'name="AVW_Para_PM" surface_pair="AVW_Para_PM"')
         .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>')
         .replace('<two_pass>1</two_pass>', f'<two_pass>{tp}</two_pass>')
         .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>')
         .replace('<offset>2.0</offset>', '<offset>0</offset>'))
    for s in ('AVW_Para_PM', f'<penalty>{pen:g}</penalty>', f'<two_pass>{tp}</two_pass>', f'<search_radius>{sr:g}</search_radius>',
              '<offset>0</offset>'):
        assert s in c, s
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE: contact AVW_Para_PM (sliding-elastic PENALTY {pen:g} auto, two_pass {tp}, search_radius {sr:g}, '
           f'offset 0, seg_up 2, shell_bottom 0/0): primary = {len(prim)} facets of {" + ".join(SHEETS)} wound to face the PM '
           f'({flipped} rewound), secondary = {len(pm)} PM_Plane facets reversed (facing the sheets). Sheet facet centres vs '
           f'nearest PM facet at rest: {(side < 0).sum()} on the -normal side, {(side > 0).sum()} on the + side '
           f'(closest {side[side < 0].max():+.2f} mm)')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch296.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--sr', 1.0), g('--pen', 0.5), int(g('--tp', 0)))
