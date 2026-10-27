"""Batch C023 (2026-10-26, Stephen's long session). curve_cross_diag.py + tri_check on C017: the remaining deep curve
crossings (0.5-0.85 mm, mid-canal, rest z ~5-8 and the strap corner) are curve-skin edges (curve_all) passing through faces
whose nodes are all in walls_LA_primary, i.e. the wall's outer surface where the curve ends are tied into the walls - a
curve-vs-WALL crossing. The model has no curve-vs-wall contact (C018 / C022 aimed at the curve end faces and did nothing).
Here: a new contact curve_wall, primary = the curve skin (curve_all), secondary = the faces of walls_LA_primary whose
centroid lies within RADIUS mm (rest) of a curve_all node (new surface walls_near_curve), sliding-elastic with the
curve self-contact's settings (penalty 5, search radius 0.5, one pass). NOT IN SOURCE.
usage: py -3.10 build_batchC023.py SRC_RUN NEW_RUN [--radius 3]"""
import argparse
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from paths import RUNS_DIR  # noqa: E402


def convert(t, radius=3.0):
    X = {int(i): np.array([float(v) for v in xyz.split(',')])
         for i, xyz in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    cb = re.search(r'<Surface name="curve_all">(.*?)</Surface>', t, re.S).group(1)
    cn = sorted({int(v) for f in re.findall(r'>([^<]+)</quad4>', cb) for v in f.split(',')})
    tree = cKDTree(np.array([X[n] for n in cn]))
    wb = re.search(r'<Surface name="walls_LA_primary">(.*?)</Surface>', t, re.S).group(1)
    faces = re.findall(r'<(quad4|tri3) id="\d+">([^<]+)</\1>', wb)
    keep = [(k, f) for k, f in faces
            if tree.query(np.mean([X[int(v)] for v in f.split(',')], 0))[0] < radius]
    assert keep, 'no wall faces near the curve'
    surf = ('\t\t<Surface name="walls_near_curve">\n'
            + ''.join(f'\t\t\t<{k} id="{i + 1}">{f}</{k}>\n' for i, (k, f) in enumerate(keep)) + '\t\t</Surface>\n')
    i = t.index('\t\t<Surface name="curve_all">')
    t = t[:i] + surf + t[i:]
    pair = ('\t\t<SurfacePair name="curve_wall">\n\t\t\t<primary>curve_all</primary>\n'
            '\t\t\t<secondary>walls_near_curve</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name="curve_all_self">')
    t = t[:i] + pair + t[i:]
    c = re.search(r'\t\t<contact name="curve_all_self" surface_pair="curve_all_self" type="sliding-elastic">.*?</contact>\n',
                  t, re.S).group(0)
    t = t.replace(c, c + c.replace('name="curve_all_self" surface_pair="curve_all_self"',
                                   'name="curve_wall" surface_pair="curve_wall"'))
    return t, [f'contact curve_wall: curve skin (curve_all, {len(cn)} nodes) against the {len(keep)} of {len(faces)} '
               f'walls_LA_primary faces within {radius:g} mm of the curve (new surface walls_near_curve), settings as '
               f'curve_all_self (penalty 5, search radius 0.5, one pass) [NOT IN SOURCE]']


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('new'); ap.add_argument('--radius', type=float, default=3)
    a = ap.parse_args()
    t = open(os.path.join(RUNS_DIR, a.src, a.src + '.feb'), encoding='utf-8', newline='').read()
    t, log = convert(t, a.radius)
    print(a.new); print('\n'.join(log))
    write(a.new, t, f'base: claude_diag/runs/{a.src}/{a.src}.feb (tools/build_batchC023.py)\nONE change (the remaining mid-canal '
                    'crossing is curve skin through the wall surface where the curve ends are tied in):\n- '
                    + '\n- '.join(log) + '\n')
