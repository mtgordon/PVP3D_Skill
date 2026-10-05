"""Resultant of each pressure load: F = -p sum(n A) over the load's facets (FEBio: a positive pressure acts against the
facet normal given by the node order). Prints the unit direction and the magnitude at the model's pressure values.
usage: py -3.10 load_direction.py MODEL.feb"""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from febmodel import Feb

f = Feb(sys.argv[1])
root = ET.parse(sys.argv[1]).getroot()
for sl in root.find('Loads'):
    if sl.get('type') != 'pressure':
        continue
    p = float(sl.find('pressure').text)
    nA = np.zeros(3)
    area = 0.0
    for t in f.surface_tris(sl.get('surface')):
        a, b, c = (np.array(f.nodes[i]) for i in t)
        v = 0.5 * np.cross(b - a, c - a)
        nA += v
        area += np.linalg.norm(v)
    F = -p * nA
    d = F / np.linalg.norm(F)
    print(f'{sl.get("name"):13s} p {p:g} MPa, area {area:6.0f} mm2: |F| {np.linalg.norm(F):6.1f} N of p*A {p * area:6.1f} N; '
          f'direction ({d[0]:+.2f}, {d[1]:+.2f}, {d[2]:+.2f})')
