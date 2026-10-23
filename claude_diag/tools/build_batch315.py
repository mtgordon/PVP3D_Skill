"""Batch 315 (2026-10-17 evening; L296j: with two_pass the PM arch's FREE inner-edge nodes chatter on the loft sheet (crawl at
t 0.301) and still 3 PM nodes go through the sheets; L296h let them through from t 0.30. Sliding-elastic searches along the
primary surface's normal, so a shell edge moving sideways into a sheet is never seen (skill gotcha 32); a solid's edge is a
real face with its own normal).
NOT IN SOURCE: the PM arch (PM_Plane, 70 quad4 shells, thickness 2 mm, PM_Yeoh) becomes a solid: NL layers of hex8 through the
2 mm, extruded from each quad toward +n (away from the AVW-Para lofts: 12 loft nodes lie 0.5-2 mm on the -n side at rest, and
the shell contact acted at the node plane; a -n extrusion overlapped the lofts and inverted elements in the first step; FEBio 4
shell nodes are the shell's TOP face with the thickness on -n, so this moves the PM's material 2 mm), so the shell nodes stay as
the solid's -n face (every linear constraint / spring anchor tied onto them unchanged). The new nodes on the
clamped outer arc join a new zero-displacement set; the shell BC on the arc ("zero shell displacement") is removed.
The contact AVW_Para_solid_PM (from build_batch312 --contact) gets the PM solid's outward faces (top, bottom, the free inner
edge and the side edges; not the clamped arc's faces) as its secondary. With --rev a second contact PM_solid_AVW_Para (PM
outward faces primary, the solid lofts' outward faces secondary; same settings) checks the PM's nodes against the lofts.
--revedge (L315c): the second contact's primary = only the PM's EDGE faces (the free inner edge and the arch's ends), not its
two big faces (L315b crawled: loft nodes bounced between the two opposing contacts on the big faces).
usage: py -3.10 build_batch315.py SRC_RUN NEW_RUN [--layers 2] [--rev | --revedge]"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

TH = 2.0


def main(src, new, nl, rev, revedge=False):
    t, X, E = load(src)
    assert '<ShellDomain name="PM_Plane" mat="PM_Yeoh">' in t and f'<shell_thickness>{TH:g}</shell_thickness>' in t
    quads = [c for nm, c in E.values() if nm == 'PM_Plane']
    arc = {int(v) for v in re.findall(r'\d+', re.search(r'<NodeSet name="BC-PM-outer-arc">(.*?)</NodeSet>', t, re.S).group(1))}
    nsum = {}
    for q in quads:
        a, b, c, d = q
        n = np.cross(X[c] - X[a], X[d] - X[b])
        for v in q:
            nsum[v] = nsum.get(v, 0) + n
    nid = max(X) + 1; eid = max(E) + 1
    lay = {0: {v: v for v in nsum}}; rows = []
    for k in range(1, nl + 1):
        lay[k] = {}
        for v, n in sorted(nsum.items()):
            p = X[v] + (TH * k / nl) * n / np.linalg.norm(n)
            lay[k][v] = nid; rows.append(f'\t\t\t<node id="{nid}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>\n'); nid += 1
    P = {v: X[v] for v in nsum}
    for k in range(1, nl + 1):
        for v, n in nsum.items():
            P[lay[k][v]] = X[v] + (TH * k / nl) * n / np.linalg.norm(n)
    hexes = []
    for q in quads:
        for k in range(nl):
            lo, hi = lay[k], lay[k + 1]
            h = [lo[v] for v in q] + [hi[v] for v in q]
            hexes.append(h)
    # volume check (5-tet split)
    T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]
    vmin = min(np.dot(np.cross(P[h[b]] - P[h[a]], P[h[c]] - P[h[a]]), P[h[d]] - P[h[a]]) for h in hexes for a, b, c, d in T5)
    assert vmin > 0, f'inverted hex at rest (min tet volume {vmin})'
    # remove the shells
    for pat in (r'\t\t<Elements type="quad4" name="PM_Plane">.*?</Elements>\n', r'\t\t<ShellDomain name="PM_Plane".*?</ShellDomain>\n',
                r'\t\t<bc name="BC-PM-outer-arc-shell".*?</bc>\n'):
        t, k = re.subn(pat, '', t, count=1, flags=re.S); assert k == 1, pat
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="PM_solid_layers">\n' + ''.join(rows) + '\t\t</Nodes>\n' + t[i:]
    i = t.rindex('</Elements>') + len('</Elements>\n')
    t = t[:i] + '\t\t<Elements type="hex8" name="PM_Plane">\n' + ''.join(
        f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hexes)) + '\t\t</Elements>\n' + t[i:]
    i = t.index('\n', t.rindex('<SolidDomain ')) + 1
    t = t[:i] + '\t\t<SolidDomain name="PM_Plane" mat="PM_Yeoh" />\n' + t[i:]
    newarc = sorted(lay[k][v] for k in range(1, nl + 1) for v in arc)
    i = t.index('<NodeSet name="BC-PM-outer-arc">')
    t = t[:i] + '<NodeSet name="BC-PM-outer-arc-solid">' + ','.join(map(str, newarc)) + '</NodeSet>\n\t\t' + t[i:]
    b = re.search(r'\t\t<bc name="BC-PM-outer-arc" node_set="BC-PM-outer-arc" type="zero displacement">.*?</bc>\n', t, re.S).group(0)
    t = t.replace(b, b + b.replace('BC-PM-outer-arc', 'BC-PM-outer-arc-solid'), 1)
    # outward faces of the PM solid, without the faces on the clamped arc
    FI = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    cnt = Counter(frozenset(h[j] for j in f) for h in hexes for f in FI)
    arcall = arc | set(newarc)
    faces, side = [], []
    for h in hexes:
        for k, f in enumerate(FI):
            q = [h[j] for j in f]
            if cnt[frozenset(q)] == 1 and not set(q) <= arcall:
                faces.append(q)
                if k >= 2:
                    side.append(q)   # the PM's edge faces (not its two big faces)
    fx = ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, f))}</quad4>\n' for k, f in enumerate(faces, 1))
    m = re.search(r'(\t\t<Surface name="AVW_Para_solid_PM_secondary">\n)(.*?)(\t\t</Surface>\n)', t, re.S)
    assert m, 'source needs build_batch312 --contact'
    t = t.replace(m.group(0), m.group(1) + fx + m.group(3), 1)
    msg = (f'NOT IN SOURCE: PM arch as a solid: {len(hexes)} hex8 ({nl} layers through the {TH:g} mm, extruded toward +n, away from the lofts; the shell '
           f'nodes stay as the -n face), PM_Yeoh; {len(newarc)} new outer-arc nodes in BC-PM-outer-arc-solid (zero displacement); '
           f'BC-PM-outer-arc-shell removed; AVW_Para_solid_PM secondary = the PM solid\'s {len(faces)} outward faces (clamped-arc faces left out)')
    if rev:
        prim = re.search(r'\t\t<Surface name="AVW_Para_solid_PM_primary">\n(.*?)\t\t</Surface>\n', t, re.S).group(1)
        fxe = ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, f))}</quad4>\n' for k, f in enumerate(side, 1)) if revedge else fx
        s = ('\t\t<Surface name="PM_solid_AVW_Para_primary">\n' + fxe + '\t\t</Surface>\n\t\t<Surface name="PM_solid_AVW_Para_secondary">\n'
             + prim + '\t\t</Surface>\n\t\t<SurfacePair name="PM_solid_AVW_Para">\n\t\t\t<primary>PM_solid_AVW_Para_primary</primary>\n'
             '\t\t\t<secondary>PM_solid_AVW_Para_secondary</secondary>\n\t\t</SurfacePair>\n')
        i = t.index('\t\t<SurfacePair name=')
        t = t[:i] + s + t[i:]
        c = re.search(r'\t\t<contact name="AVW_Para_solid_PM".*?</contact>\n', t, re.S).group(0)
        t = t.replace(c, c + c.replace('name="AVW_Para_solid_PM" surface_pair="AVW_Para_solid_PM"', 'name="PM_solid_AVW_Para" surface_pair="PM_solid_AVW_Para"'), 1)
        msg += f'; + contact PM_solid_AVW_Para (PM {"EDGE faces only (" + str(len(side)) + ")" if revedge else "outward faces"} primary vs the solid lofts\' outward faces, settings as AVW_Para_solid_PM)'
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch315.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], int(a[a.index('--layers') + 1]) if '--layers' in a else 2, '--rev' in a or '--revedge' in a, '--revedge' in a)
