"""Batch 189 (2026-10-06, the user: "yes, build it and start it"): the tube's curves carry no pressure.
In every tube run the intraabdominal pressure (0.014 MPa, lc 1) acts only on the source's surfaces (Load-AVW / Load-PVW /
Load-PeB-top / Load-top, + Load-LA); the curves (canal_wrap_hex_avw / _pvw; NOT IN SOURCE geometry) are pushed only through
the 536 nodes they share with the walls. One change on the working unrefined tube run L169_wrapfull_endspr_selfall_walls035_P1:
  L189_selfall_curvepress_P1   + a pressure load "Load-curves" (0.014 MPa, lc 1, as Load-AVW) on the curves' outer faces that
                               face the same way as the loaded wall next to them (NOT IN SOURCE)
Face choice: each boundary face of the curves (surface curve_all), oriented outward from its own hex; kept if the nearest
Load-AVW / Load-PVW face (by centroid) is within MAXD mm and the two outward normals agree (cos > COS). Written with the same
winding sense as the Load-AVW faces (checked against their own hexes).
usage: py -3.10 build_batch189.py [SRC NAME]   (default L169 -> L189_selfall_curvepress_P1; 2026-10-06: also
       L187_selfall_botfull_P1 -> L190_botfull_curvepress_P1, the user's equal end attachments + the curve pressure)"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L169_wrapfull_endspr_selfall_walls035_P1'
NAME = 'L189_selfall_curvepress_P1'
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
MAXD, COS = 6.0, 0.5


def faces_of(surf):
    return [[int(v) for v in f.text.split(',')] for f in surf]


if __name__ == '__main__':
    if len(sys.argv) == 3:
        SRC, NAME = sys.argv[1], sys.argv[2]
    p = os.path.join(RUNS, SRC, SRC + '.feb')
    t = open(p, encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    owner = {}          # node -> list of (domain, hex centroid)
    for e in mesh.findall('Elements'):
        for el in e:
            ids = [int(v) for v in el.text.split(',')]
            c = np.mean([X[i] for i in ids], axis=0)
            for i in ids:
                owner.setdefault(i, []).append((e.get('name'), frozenset(ids), c))
    surfs = {s.get('name'): faces_of(s) for s in mesh.findall('Surface')}

    def hexc(f, doms):
        cand = [(fs, c) for (d, fs, c) in owner[f[0]] if d in doms and all(v in fs for v in f)]
        assert cand, f
        return cand[0][1]

    def geom(f, doms):
        P = np.array([X[v] for v in f])
        c = P.mean(0)
        n = np.cross(P[1] - P[0], P[-1] - P[0])
        n /= np.linalg.norm(n)
        out = c - hexc(f, doms)
        return c, n, (1 if n @ out > 0 else -1)

    walls = ('_PickedSet347', '_PickedSet64', '_PickedSet346')
    L = [geom(f, walls) for f in surfs['Load-AVW_surf'] + surfs['Load-PVW_surf']]
    wsense = {s for _, _, s in L}
    assert len(wsense) == 1, f'Load-AVW / PVW faces wound both ways: {wsense}'
    wsense = wsense.pop()
    LC = np.array([c for c, _, _ in L])
    LN = np.array([n * s for _, n, s in L])        # outward normals
    keep = []
    for f in surfs['curve_all']:
        c, n, s = geom(f, CURVES)
        nout = n * s
        d = np.linalg.norm(LC - c, axis=1)
        j = int(d.argmin())
        if d[j] <= MAXD and nout @ LN[j] > COS:
            keep.append(f if s == wsense else [f[0]] + f[1:][::-1])
    assert keep
    tag = lambda f: 'quad4' if len(f) == 4 else 'tri3'
    surf = '\t\t<Surface name="Load-curves_surf">\n' + ''.join(
        f'\t\t\t<{tag(f)} id="{i}">{",".join(map(str, f))}</{tag(f)}>\n' for i, f in enumerate(keep, 1)) + '\t\t</Surface>\n'
    i = t.index('\t\t<Surface name="curve_all">')
    t = t[:i] + surf + t[i:]
    a = t.index('<surface_load name="Load-AVW"')
    b = t.index('</surface_load>', a) + len('</surface_load>')
    blk = t[a:b].replace('name="Load-AVW"', 'name="Load-curves"').replace('surface="Load-AVW_surf"', 'surface="Load-curves_surf"')
    t = t[:b] + '\n\t\t' + blk + t[b:]
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch189.py)\n'
                   f'NOT IN SOURCE: pressure "Load-curves" (as Load-AVW: 0.014 MPa, lc 1) on {len(keep)} of '
                   f'{len(surfs["curve_all"])} curve boundary faces (outer faces within {MAXD} mm of a Load-AVW / PVW face and '
                   f'facing the same way, cos > {COS})\n')
    print(len(keep), 'of', len(surfs['curve_all']), 'curve faces loaded')
