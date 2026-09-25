"""Lightweight FEBio .feb parser for structural diagnostics (no FEBio needed)."""
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter

import numpy as np


class Feb:
    def __init__(self, path):
        self.path = path
        self.tree = ET.parse(path)
        self.root = self.tree.getroot()
        self.nodes = {}            # id -> np.array(3)
        self.node_block = {}       # id -> Nodes block name
        self.elem_blocks = {}      # name -> (type, {eid: [nids]})
        self.elem_owner = {}       # eid -> block name
        self.nodesets = {}         # name -> [nids]
        self.surfaces = {}         # name -> [(type, [nids])]
        self.surfpairs = {}        # name -> (primary, secondary)
        self.discsets = {}         # name -> [(n1, n2)]
        self._parse_mesh()

    def _parse_mesh(self):
        mesh = self.root.find('Mesh')
        for child in mesh:
            tag = child.tag
            name = child.get('name')
            if tag == 'Nodes':
                for n in child:
                    nid = int(n.get('id'))
                    self.nodes[nid] = np.array([float(v) for v in n.text.split(',')])
                    self.node_block[nid] = name
            elif tag == 'Elements':
                et = child.get('type')
                d = {}
                for e in child:
                    eid = int(e.get('id'))
                    d[eid] = [int(v) for v in e.text.split(',')]
                    self.elem_owner[eid] = name
                self.elem_blocks[name] = (et, d)
            elif tag == 'NodeSet':
                txt = (child.text or '').strip()
                ids = []
                if txt:
                    ids = [int(v) for v in re.split(r'[,\s]+', txt) if v]
                for n in child.findall('n'):
                    ids.append(int(n.get('id')))
                self.nodesets[name] = ids
            elif tag == 'Surface':
                facets = []
                for f in child:
                    facets.append((f.tag, [int(v) for v in f.text.split(',')]))
                self.surfaces[name] = facets
            elif tag == 'SurfacePair':
                self.surfpairs[name] = (child.find('primary').text.strip(),
                                        child.find('secondary').text.strip())
            elif tag == 'DiscreteSet':
                self.discsets[name] = [tuple(int(v) for v in d.text.split(','))
                                       for d in child.findall('delem')]

    # --- helpers -------------------------------------------------------
    def domain_nodes(self, name):
        return sorted({n for conn in self.elem_blocks[name][1].values() for n in conn})

    def nodes_in_elements(self):
        s = set()
        for et, d in self.elem_blocks.values():
            for conn in d.values():
                s.update(conn)
        return s

    def surface_nodes(self, name):
        return sorted({n for _, conn in self.surfaces[name] for n in conn})

    def surface_tris(self, name):
        """Split each facet into triangles (quad4 -> 2 tris)."""
        tris = []
        for t, c in self.surfaces[name]:
            if len(c) == 3:
                tris.append(c)
            elif len(c) == 4:
                tris.append([c[0], c[1], c[2]])
                tris.append([c[0], c[2], c[3]])
            else:
                tris.append(c[:3])
        return tris

    def section(self, tag):
        return self.root.find(tag)


def pt_tri_dist(p, a, b, c):
    """Closest distance from point p to triangle abc (Ericson)."""
    ab = b - a; ac = c - a; ap = p - a
    d1 = ab @ ap; d2 = ac @ ap
    if d1 <= 0 and d2 <= 0:
        return np.linalg.norm(p - a)
    bp = p - b
    d3 = ab @ bp; d4 = ac @ bp
    if d3 >= 0 and d4 <= d3:
        return np.linalg.norm(p - b)
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return np.linalg.norm(p - (a + v * ab))
    cp = p - c
    d5 = ab @ cp; d6 = ac @ cp
    if d6 >= 0 and d5 <= d6:
        return np.linalg.norm(p - c)
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return np.linalg.norm(p - (a + w * ac))
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return np.linalg.norm(p - (b + w * (c - b)))
    denom = 1.0 / (va + vb + vc)
    v = vb * denom; w = vc * denom
    return np.linalg.norm(p - (a + ab * v + ac * w))


def dist_to_surface(feb, pts, tris):
    """Min distance from each point to a triangle soup; returns (dist, tri_index)."""
    A = np.array([feb.nodes[t[0]] for t in tris])
    B = np.array([feb.nodes[t[1]] for t in tris])
    C = np.array([feb.nodes[t[2]] for t in tris])
    cen = (A + B + C) / 3
    out = []
    for p in pts:
        # prefilter by centroid distance
        dc = np.linalg.norm(cen - p, axis=1)
        idx = np.argsort(dc)[:60]
        best = (1e30, -1)
        for i in idx:
            d = pt_tri_dist(p, A[i], B[i], C[i])
            if d < best[0]:
                best = (d, i)
        out.append(best)
    return out


if __name__ == '__main__':
    f = Feb(sys.argv[1])
    print(len(f.nodes), 'nodes;', len(f.elem_blocks), 'element blocks')
