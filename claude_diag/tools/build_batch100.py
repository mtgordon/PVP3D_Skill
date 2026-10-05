"""Batch 100 (2026-09-28, the user away): directional or local wall stiffness instead of stiffer walls everywhere (the
user's follow-up to batch 99: "stiffness in one direction only ... e.g. fibres across the canal in the AVW or a stiffer band
along the seam"). NOT IN SOURCE: the source's walls are isotropic (Marlow). Syntax checked on
claude_diag/fibre_mini/fibre_mini.py (a fibre adds stiffness only along its direction; per-element directions work through a
mat_axis ElementData with the fibre vector 1,0,0 in the element's axes, but a fibre without its <fiber> crashes silently;
a mapped c1 works; a mixture's mass comes from its own top-level density, not its solids').

fibres(ef, walls): each wall material -> an uncoupled solid mixture of its faithful Yeoh (same c1, c2, k, density) + one
    fiber-pow-linear-uncoupled family (E = ef, beta 2, lam0 1.01: about linear with modulus ef from the start, tension only)
    running across the canal: the global x axis (lateral), which lies in the wall's tangent plane in every AVW element
    (checked per element against the nearest canal-contact facet's normal: within 0.1 deg; the canal surfaces are flat
    across x); the fibre adds ~ef to the small-strain modulus along that direction only. ef =
    6 c1 (f - 1) matches the c1 x f law along the fibre (ef 0.0358 = the c1 x2 law's extra small-strain stiffness).
band(width, factor, walls): the Yeoh c1 of every element whose centroid lies within `width` mm of a lateral canal-seam
    node (build_batch90.seam_nodes: the 116 AVW / cervix edge nodes) x factor, through a mapped c1 (the rest faithful).
usage: py -3.10 build_batch100.py NAME ...   (see BUILDS; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from variants7 import Model7, JOBS, emit
from build_batch90 import seam_nodes

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
DOMS = {'Vagina_AVW': '_PickedSet347', 'Vagina_PVW': '_PickedSet64', 'Vagina_Cervix': '_PickedSet346'}
C1 = 0.005964127


def _mesh_data(m):
    md = m.root.find('MeshData')
    if md is None:
        md = ET.Element('MeshData')
        kids = list(m.root)
        m.root.insert(kids.index(m.root.find('MeshDomains')) + 1, md)
    return md


def _elements(m, dom):
    X = m.nodes()
    blk = m.elem_blocks()[dom]
    ids = [int(e.get('id')) for e in blk]
    cen = np.array([np.mean([X[int(v)] for v in e.text.split(',')], axis=0) for e in blk])
    return ids, cen


def _facets(m, name):
    X = m.nodes()
    s = next(s for s in m.mesh.findall('Surface') if s.get('name') == name)
    F = [[X[int(v)] for v in f.text.split(',')] for f in s]
    cen = np.array([np.mean(f, axis=0) for f in F])
    nrm = np.array([np.cross(f[2] - f[0], f[-1] - f[1]) if len(f) == 4 else np.cross(f[1] - f[0], f[2] - f[0]) for f in F])
    return cen, nrm / np.linalg.norm(nrm, axis=1)[:, None]


def across_dirs(m, dom):
    """per element: x projected onto the wall's tangent plane (normal of the nearest canal facet of either wall)."""
    ids, cen = _elements(m, dom)
    fc, fn = map(np.vstack, zip(_facets(m, 'SlidingElastic1Secondary'), _facets(m, 'SlidingElastic1Primary')))
    ex = np.array([1.0, 0, 0])
    A, N, fall = [], [], 0
    for c in cen:
        n = fn[np.argmin(np.linalg.norm(fc - c, axis=1))]
        a = ex - (ex @ n) * n
        if np.linalg.norm(a) < 0.5:                  # the wall faces sideways: keep x (the tangent plane has no x)
            fall += 1
            a = ex
        A.append(a / np.linalg.norm(a))
        N.append(n if abs(n @ A[-1]) < 0.99 else np.array([0, 0, 1.0]))
    return ids, np.array(A), np.array(N), fall


def fibres(ef, walls=('Vagina_AVW',)):
    def fn(m):
        for w in walls:
            mat = m._material(w)
            assert mat.get('type') == 'Yeoh', w
            p = {c.tag: c.text for c in mat}
            for c in list(mat):
                mat.remove(c)
            mat.set('type', 'uncoupled solid mixture')
            ET.SubElement(mat, 'density').text = p['density']
            ET.SubElement(mat, 'k').text = p['k']
            s = ET.SubElement(mat, 'solid', type='Yeoh')
            for k in ('c1', 'c2'):
                ET.SubElement(s, k).text = p[k]
            s = ET.SubElement(mat, 'solid', type='fiber-pow-linear-uncoupled')
            ET.SubElement(s, 'fiber', type='vector').text = '1,0,0'
            ET.SubElement(s, 'E').text = '%.6g' % ef
            ET.SubElement(s, 'beta').text = '2'
            ET.SubElement(s, 'lam0').text = '1.01'
            ids, A, N, fall = across_dirs(m, DOMS[w])
            dev = float(np.degrees(np.arccos(np.clip(np.abs(A[:, 0]).min(), 0, 1))))
            assert dev < 1, f'{w}: the across-canal direction departs {dev:.1f} deg from x: write a mat_axis ElementData'
            m.log.append(f'NOT IN SOURCE (directional wall stiffness): {w} = its faithful Yeoh (c1 {p["c1"]}, c2 {p["c2"]}, '
                         f'k {p["k"]}) + one fiber-pow-linear-uncoupled family (E {ef:.4g} MPa, beta 2, lam0 1.01; tension '
                         f'only) along x = across the canal (x lies in the wall\'s tangent plane in all {len(ids)} elements, '
                         f'within {dev:.2f} deg: the canal surfaces are flat across x)')
    return fn


def fibres_along(ef, walls=('Vagina_AVW',)):
    """as fibres(), but the fibre runs ALONG the canal: per element the wall's tangent perpendicular to x (n x ex, n =
    the nearest canal facet's normal), written as a mat_axis ElementData (a = that direction, d = n) with the fibre vector
    1,0,0 in the element's axes (checked on claude_diag/fibre_mini/fibre_mini.py)."""
    def fn(m):
        md = _mesh_data(m)
        for w in walls:
            mat = m._material(w)
            assert mat.get('type') == 'Yeoh', w
            p = {c.tag: c.text for c in mat}
            for c in list(mat):
                mat.remove(c)
            mat.set('type', 'uncoupled solid mixture')
            ET.SubElement(mat, 'density').text = p['density']
            ET.SubElement(mat, 'k').text = p['k']
            s = ET.SubElement(mat, 'solid', type='Yeoh')
            for k in ('c1', 'c2'):
                ET.SubElement(s, k).text = p[k]
            s = ET.SubElement(mat, 'solid', type='fiber-pow-linear-uncoupled')
            ET.SubElement(s, 'fiber', type='vector').text = '1,0,0'
            ET.SubElement(s, 'E').text = '%.6g' % ef
            ET.SubElement(s, 'beta').text = '2'
            ET.SubElement(s, 'lam0').text = '1.01'
            ids, A, N, fall = across_dirs(m, DOMS[w])
            ed = ET.SubElement(md, 'ElementData', type='mat_axis', elem_set=DOMS[w])
            ang = []
            for i, n in enumerate(N, 1):
                a = np.cross(n, [1.0, 0, 0])
                a /= np.linalg.norm(a)
                ang.append(np.degrees(np.arccos(abs(a[2]))))
                e = ET.SubElement(ed, 'elem', lid=str(i))
                ET.SubElement(e, 'a').text = ','.join('%.6f' % v for v in a)
                ET.SubElement(e, 'd').text = ','.join('%.6f' % v for v in n)
            m.log.append(f'NOT IN SOURCE (directional wall stiffness): {w} = its faithful Yeoh (c1 {p["c1"]}, c2 {p["c2"]}, '
                         f'k {p["k"]}) + one fiber-pow-linear-uncoupled family (E {ef:.4g} MPa, beta 2, lam0 1.01; tension '
                         f'only) ALONG the canal: per element the wall tangent perpendicular to x (n x ex, n = the nearest '
                         f'canal facet normal; angle to z median {np.median(ang):.0f} deg, range {min(ang):.0f}-'
                         f'{max(ang):.0f}), a mat_axis ElementData on {len(ids)} elements')
    return fn


def band(width, factor, walls=('Vagina_AVW', 'Vagina_Cervix', 'Vagina_PVW')):
    def fn(m):
        seam, _, X = seam_nodes(m)
        S = np.array([X[v] for v in seam])
        md = _mesh_data(m)
        for w in walls:
            mat = m._material(w)
            assert mat.get('type') == 'Yeoh', w
            c1 = float(mat.find('c1').text)
            ids, cen = _elements(m, DOMS[w])
            d = np.min(np.linalg.norm(cen[:, None, :] - S[None, :, :], axis=2), axis=1)
            inb = d < width
            name = f'{w}_c1map'
            mat.find('c1').set('type', 'map')
            mat.find('c1').text = name
            ed = ET.SubElement(md, 'ElementData', name=name, elem_set=DOMS[w])
            for i, b in enumerate(inb, 1):
                ET.SubElement(ed, 'elem', lid=str(i)).text = '%.7g' % (c1 * (factor if b else 1))
            m.log.append(f'NOT IN SOURCE (a stiffer band along the canal seam): {w} Yeoh c1 x{factor} ({c1 * factor:.7g}) in '
                         f'the {int(inb.sum())} of {len(ids)} elements whose centroid lies within {width} mm of a lateral '
                         f'canal-seam node ({len(seam)} nodes); the rest keep c1 {c1} (a mapped c1)')
    return fn


BASE = 'L85_newline_vwyeoh_rhoi0'
BUILDS = (('L100_newline_vwyeoh_rhoi0_avwfib2', BASE, fibres(6 * C1 * (2 - 1))),
          ('L100_newline_vwyeoh_rhoi0_avwfib3', BASE, fibres(6 * C1 * (3 - 1))),
          ('L100_newline_vwyeoh_rhoi0_band5x3', BASE, band(5, 3)),
          # 16:40: the across-canal fibres crawl at t 0.75-0.79 like the faithful walls, while the AVW alone at c1 x2
          # (isotropic) passed 0.95: the stiffness that matters is not across the canal. Along the canal instead:
          ('L102_newline_vwyeoh_rhoi0_avwfibalong2', BASE, fibres_along(6 * C1 * (2 - 1))),
          # 17:07: the AVW alone at isotropic c1 x2 crawled at t 0.96 (canal seam) where all three walls at x2 reached
          # t = 1: the PVW / cervix matter near full load. The along-canal fibres in all three walls:
          ('L102_newline_vwyeoh_rhoi0_allfibalong2', BASE,
           fibres_along(6 * C1 * (2 - 1), ('Vagina_AVW', 'Vagina_PVW', 'Vagina_Cervix'))),
          # 17:58: the 5 mm band (c1 x3 in a third of the wall elements) reached t = 1 with Ba at t 0.76 like c1 x1.5
          # everywhere. A narrower band (3 mm: AVW 208, cervix 103, PVW 198 elements, ~10 %):
          ('L102_newline_vwyeoh_rhoi0_band3x3', BASE, band(3, 3)),
          # 18:02: and the 5 mm band at x2 (does a softer band suffice?):
          ('L102_newline_vwyeoh_rhoi0_band5x2', BASE, band(5, 2)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
