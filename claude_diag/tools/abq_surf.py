"""Abaqus element-based surfaces -> facet coordinates, and FEBio surface matching by coordinate.

part_elements(part)          {eid: [local node ids]} for every *Element block of a part
elset(name)                  element ids of an (assembly or part) *Elset, 'generate' aware
surface_faces(surf)          [(elset, face)] rows of an assembly *Surface, type=ELEMENT
face_nodes(conn, face)       node ids of an Abaqus face label (C3D8 S1..S6, shell SPOS/SNEG)
"""
import re

import numpy as np

from inpmap import INP, part_nodes

# Abaqus C3D8 face -> local node indices (outward normal by right-hand rule)
HEX_FACES = {'S1': (0, 1, 2, 3), 'S2': (4, 7, 6, 5), 'S3': (0, 4, 5, 1),
             'S4': (1, 5, 6, 2), 'S5': (2, 6, 7, 3), 'S6': (3, 7, 4, 0)}

_LINES = None


def lines():
    global _LINES
    if _LINES is None:
        _LINES = open(INP, encoding='latin-1').read().splitlines()
    return _LINES


def part_elements(part):
    out = {}
    inpart = inel = False
    for ln in lines():
        s = ln.strip()
        if s.lower().startswith('*part,'):
            m = re.search(r'name=("?)(.+?)\1\s*$', s)
            inpart = bool(m) and m.group(2) == part
            continue
        if not inpart:
            continue
        if s.lower().startswith('*end part'):
            break
        if s.startswith('*'):
            inel = s.lower().startswith('*element,')
            continue
        if inel and s:
            v = [int(x) for x in s.replace(' ', '').split(',') if x]
            out[v[0]] = v[1:]
    return out


def elset(name, instance=None):
    """Element ids of *Elset elset=name (first match; optionally require instance=)."""
    out = []
    grab = gen = False
    for ln in lines():
        s = ln.strip()
        if s.startswith('*'):
            if grab:
                break
            low = s.lower()
            if low.startswith('*elset') and re.search(r'elset=%s\s*(,|$)' % re.escape(name), s):
                if instance and ('instance=%s' % instance) not in s:
                    continue
                grab = True
                gen = 'generate' in low
            continue
        if grab and s:
            v = [int(x) for x in s.replace(' ', '').split(',') if x]
            out += list(range(v[0], v[1] + 1, v[2] if len(v) > 2 else 1)) if gen else v
    return out


def surface_faces(surf):
    rows = []
    grab = False
    for ln in lines():
        s = ln.strip()
        if s.startswith('*'):
            if grab:
                break
            if s.lower().startswith('*surface') and re.search(r'name=%s\s*(,|$)' % re.escape(surf), s):
                grab = True
            continue
        if grab and s:
            a, b = [t.strip() for t in s.split(',')[:2]]
            rows.append((a, b))
    return rows


def face_nodes(conn, face):
    if face in HEX_FACES:
        return [conn[i] for i in HEX_FACES[face]]
    if face in ('SPOS', 'SNEG'):
        return list(conn) if face == 'SPOS' else list(conn)[::-1]
    raise ValueError(face)


def surface_facets(surf, part):
    """[(elem id, face, [xyz...])] for an Abaqus element surface on one part."""
    X = part_nodes(part)
    E = part_elements(part)
    out = []
    for es, face in surface_faces(surf):
        for e in elset(es):
            ids = face_nodes(E[e], face)
            out.append((e, face, np.array([X[i] for i in ids])))
    return out


def key(P, nd=3):
    """Order-independent coordinate key of a facet."""
    return tuple(sorted(tuple(np.round(p, nd)) for p in P))


def normal(P):
    """Facet normal by right-hand rule (quad: diagonals cross product)."""
    if len(P) == 4:
        n = np.cross(P[2] - P[0], P[3] - P[1])
    else:
        n = np.cross(P[1] - P[0], P[2] - P[0])
    return n / np.linalg.norm(n)
