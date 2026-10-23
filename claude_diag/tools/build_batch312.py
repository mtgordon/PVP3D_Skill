"""Batch 312 (2026-10-17; the user: the loft hits the PM almost at a right angle. Small tests (shell_tests/gen_edge.py,
gen_edge_solid.py): sliding-elastic never catches a SHELL sheet arriving edge-on (one direction, both, two_pass, roles swapped:
it passes straight through); node-on-facet / facet-on-facet fail; a thin SOLID sheet (one hex layer, 0.49 mm) is caught).
NOT IN SOURCE: the AVW-Para loft sheets (tri3 shells, AVW-Para-L/R_sheet) become solid lofts: one layer of penta6 wedges,
thickness TH (the shells' 0.49 mm), extruded from each triangle: bottom = the existing sheet nodes (so the chain springs and
the attachments are unchanged), top = new nodes offset TH along the averaged node normal. Same material (neo-Hookean).
With --contact: AVW_Para_solid_PM, sliding-elastic between every outward face of the solid lofts (top, bottom, edges; primary)
and the PM arch (PM_Plane facets reversed, facing the lofts' side), penalty PEN auto, two_pass 0, search_radius SR, offset 0.
usage: py -3.10 build_batch312.py SRC_RUN NEW_RUN [--th 0.49] [--contact] [--pen 1] [--sr 3]"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

SHEETS = ('AVW-Para-L_sheet', 'AVW-Para-R_sheet')


def main(src, new, th, contact, pen, sr):
    t, X, E = load(src)
    nid = max(X) + 1; eid = max(E) + 1
    rows_n, blocks, doms, faces, nshell = [], [], [], [], {}
    for sh in SHEETS:
        tris = [c for nm, c in E.values() if nm == sh]
        nsum = {}
        for a, b, c in tris:
            n = np.cross(X[b] - X[a], X[c] - X[a])
            for v in (a, b, c):
                nsum[v] = nsum.get(v, 0) + n
        top = {}
        for v, n in nsum.items():
            n = n / np.linalg.norm(n); p = X[v] + th * n
            top[v] = nid; rows_n.append(f'\t\t\t<node id="{nid}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>\n'); nid += 1
        TP = {v: X[v] + th * (n / np.linalg.norm(n)) for v, n in nsum.items()}

        def ok(a, b, c):   # all six corner tets of the wedge positive (FEBio rejects a twisted wedge at initialisation)
            A, B, C, A2, B2, C2 = X[a], X[b], X[c], TP[a], TP[b], TP[c]
            vol = lambda p0, p1, p2, p3: np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0)   # noqa: E731
            tets = [(A, B, C, A2), (B, C, A, B2), (C, A, B, C2), (A2, C2, B2, A), (B2, A2, C2, B), (C2, B2, A2, C)]
            return min(vol(*q) for q in tets) > 1e-6 * th ** 3
        shells = {tuple(f) for f in tris if not ok(*f)}
        nshell[sh] = len(shells)
        tris = [f for f in tris if tuple(f) not in shells]
        pens = [(a, b, c, top[a], top[b], top[c]) for a, b, c in tris]
        blocks.append(f'\t\t<Elements type="penta6" name="{sh}_solid">\n' +
                      ''.join(f'\t\t\t<elem id="{eid + k}">{",".join(map(str, p))}</elem>\n' for k, p in enumerate(pens)) + '\t\t</Elements>\n')
        eid += len(pens)
        mat = re.search(rf'<ShellDomain name="{sh}" mat="([^"]+)"', t).group(1)
        doms.append(f'\t\t<SolidDomain name="{sh}_solid" mat="{mat}" />\n')
        # outward boundary faces
        ecnt = Counter(frozenset(e) for a, b, c in tris for e in ((a, b), (b, c), (c, a)))
        for a, b, c in tris:
            faces.append(('tri3', [a, c, b]))                     # bottom (outward = -n)
            faces.append(('tri3', [top[a], top[b], top[c]]))      # top (outward = +n)
            for u, w in ((a, b), (b, c), (c, a)):
                if ecnt[frozenset((u, w))] == 1:                  # sheet boundary edge -> side face, outward
                    q = [u, w, top[w], top[u]]
                    cen = (X[a] + X[b] + X[c]) / 3
                    nq = np.cross(X[w] - X[u], (X[u] + th * (nsum[u] / np.linalg.norm(nsum[u]))) - X[u])
                    if np.dot(nq, (X[u] + X[w]) / 2 - cen) < 0:
                        q = q[::-1]
                    faces.append(('quad4', q))
        # remove the shells, except the triangles that cannot form a valid wedge (kept as shells, same domain)
        if shells:
            m_ = re.search(rf'(\t\t<Elements type="tri3" name="{sh}">\n)(.*?)(\t\t</Elements>\n)', t, re.S)
            rows = [r for r in re.findall(r'\t\t\t<elem id="\d+">[0-9,]+</elem>\n', m_.group(2))
                    if tuple(int(v) for v in re.search(r'>([0-9,]+)<', r).group(1).split(',')) in shells]
            t = t.replace(m_.group(0), m_.group(1) + ''.join(rows) + m_.group(3), 1)
        else:
            for pat in (rf'\t\t<Elements type="tri3" name="{sh}">.*?</Elements>\n', rf'\t\t<ShellDomain name="{sh}".*?</ShellDomain>\n'):
                t, k = re.subn(pat, '', t, count=1, flags=re.S); assert k == 1, pat
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="AVW-Para_solid_loft_top">\n' + ''.join(rows_n) + '\t\t</Nodes>\n' + t[i:]
    i = t.rindex('</Elements>') + len('</Elements>\n')
    t = t[:i] + ''.join(blocks) + t[i:]
    i = t.index('\n', t.rindex('<SolidDomain ')) + 1
    t = t[:i] + ''.join(doms) + t[i:]
    msg = (f'NOT IN SOURCE: AVW-Para loft sheets as solid lofts: one layer of penta6, thickness {th:g} mm, extruded from each tri3 '
           f'(bottom = the sheet nodes, top = {len(rows_n)} new nodes along the averaged node normal), same materials; triangles that cannot form a valid wedge kept as shells: ' + str(nshell) + ')')
    if contact:
        pm = [c for nm, c in E.values() if nm == 'PM_Plane']
        surf = ('\t\t<Surface name="AVW_Para_solid_PM_primary">\n' +
                ''.join(f'\t\t\t<{g} id="{k}">{",".join(map(str, c))}</{g}>\n' for k, (g, c) in enumerate(faces, 1)) +
                '\t\t</Surface>\n\t\t<Surface name="AVW_Para_solid_PM_secondary">\n' +
                ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q[::-1]))}</quad4>\n' for k, q in enumerate(pm, 1)) +
                '\t\t</Surface>\n\t\t<SurfacePair name="AVW_Para_solid_PM">\n\t\t\t<primary>AVW_Para_solid_PM_primary</primary>\n'
                '\t\t\t<secondary>AVW_Para_solid_PM_secondary</secondary>\n\t\t</SurfacePair>\n')
        i = t.index('\t\t<SurfacePair name=')
        t = t[:i] + surf + t[i:]
        a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
        c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', 'name="AVW_Para_solid_PM" surface_pair="AVW_Para_solid_PM"')
             .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
             .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
        for s_ in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>'):
            assert s_ in c, s_
        t = t[:b] + '\t\t' + c + t[b:]
        msg += (f'; contact AVW_Para_solid_PM: {len(faces)} outward faces of the solid lofts vs PM_Plane reversed (sliding-elastic, '
                f'penalty {pen:g} auto, two_pass 0, search_radius {sr:g}, offset 0, seg_up 2)')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch312.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--th', 0.49), '--contact' in a, g('--pen', 1.0), g('--sr', 3.0))
