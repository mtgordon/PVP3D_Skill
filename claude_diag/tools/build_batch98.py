"""Batch 98 (2026-09-28 ~12:30): a canal seam that can slide but not open (the user: "Test a seam restraint that allows
sliding but not opening"; NOT IN SOURCE).

The seam springs (batch 96-97) hold the seam in every direction: they let the faithful walls converge but lift both walls,
and in the paper's cases the posterior wall stops responding to the impairment (P1 -> P2: Bp +0.6 mm with the springs,
+4.4 free, +5 in the paper). Free, the seam mainly opens (P1: up to 47 mm) rather than slides (12 mm). This restraint only
stops the opening: the canal contact's secondary facets that hold a lateral-seam node (the AVW / cervix seam strip) move into
their own sliding-elastic contact against the whole PVW / PeB surface with tension on (FEBio forum, Steve Maas: it "forces
a node from the primary surface to move onto the secondary surfaces, even when the surfaces try to separate ... displacement
along the surface is still allowed"), single pass, the canal's other settings, and an offset equal to the strip nodes'
median rest gap (without it the rest gap snaps shut in the first iteration and elements invert; checked on
claude_diag/seam_tie/seam_slide_mini.py: pulled apart, the strip stays within ~0.15 mm and slides freely).
One change each:
  L93_pen5_paperP1_seamslide              L68_pen5_paperP1 (the paper's case P1, free seam) + the sliding seam
  L93_springs_newline_rhoi0_seamslide     L87_springs_newline_rhoi0 (the new springs line + rhoi 0, Ogden walls)
  L93_newline_vwyeoh_rhoi0_seamslide      L85_newline_vwyeoh_rhoi0 (the faithful walls, free seam: crawled at t 0.78)
usage: py -3.10 build_batch98.py [NAME ...] [--gaps]   (--gaps: print the strip's rest gaps only)
"""
import copy
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from variants7 import Model7, JOBS, emit
from build_batch90 import seam_nodes

sys.path.insert(0, os.path.join(JOBS, 'claude_diag', 'seam_tie'))
from seam_mini import project  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = {a for a in sys.argv[1:] if not a.startswith('--')}


def strip_gaps(m):
    """the lateral-seam strip of the canal's secondary surface and each strip node's rest gap to the primary surface."""
    seam, fP, X = seam_nodes(m)
    surf = {s.get('name'): s for s in m.mesh.findall('Surface')}
    S = surf['SlidingElastic1Secondary']
    strip = [f for f in S if {int(v) for v in f.text.split(',')} & set(seam)]
    nodes = sorted({int(v) for f in strip for v in f.text.split(',')})
    cen = np.array([np.mean([X[n] for n in f], axis=0) for f in fP])
    gaps = []
    for v in nodes:
        near = [fP[i] for i in np.argsort(np.linalg.norm(cen - X[v], axis=1))[:12]]
        gaps.append(project(X[v], X, near)[0])
    return seam, strip, np.array(gaps)


def seam_slide(m):
    seam, strip, gaps = strip_gaps(m)
    off = float(np.median(gaps))
    mesh = m.mesh
    surf = {s.get('name'): s for s in mesh.findall('Surface')}
    S = surf['SlidingElastic1Secondary']
    n0 = len(S)
    for f in strip:
        S.remove(f)
    for i, f in enumerate(S, 1):
        f.set('id', str(i))
    new = ET.Element('Surface', {'name': 'canal_seam_strip'})
    for i, f in enumerate(strip, 1):
        q = ET.SubElement(new, f.tag, {'id': str(i)})
        q.text = f.text
    new.tail, new.text = S.tail, S.text
    kids = list(mesh)
    mesh.insert(kids.index(S) + 1, new)
    pair = next(p for p in mesh.findall('SurfacePair') if p.get('name') == 'SlidingElastic1')
    np_ = ET.Element('SurfacePair', {'name': 'canal_seam_slide'})
    ET.SubElement(np_, 'primary').text = pair.find('primary').text
    ET.SubElement(np_, 'secondary').text = 'canal_seam_strip'
    np_.tail, np_.text = pair.tail, pair.text
    mesh.insert(list(mesh).index(pair) + 1, np_)
    cont = m.root.find('Contact')
    base = next(c for c in cont if c.get('name') == 'SlidingElastic1')
    c = copy.deepcopy(base)
    c.set('name', 'canal_seam_slide')
    c.set('surface_pair', 'canal_seam_slide')
    for k, v in (('tension', 1), ('two_pass', 0), ('offset', round(off, 4))):
        c.find(k).text = str(v)
    cont.insert(list(cont).index(base) + 1, c)
    m.log.append(f'NOT IN SOURCE (a canal seam that slides but does not open): the {len(strip)} SlidingElastic1Secondary facets '
                 f'holding one of the {len(seam)} lateral-seam nodes (AVW / cervix) moved from the canal contact ({n0} -> '
                 f'{len(S)} facets) into their own sliding-elastic contact canal_seam_slide against SlidingElastic1Primary '
                 f'(PVW / PeB) with tension 1 (cannot separate, can slide), two_pass 0, offset {off:.3f} mm (the median rest '
                 f'gap of the strip\'s {len(gaps)} nodes, {gaps.min():.2f}-{gaps.max():.2f} mm), the canal\'s other settings')
    return len(strip)


def slide_set(note, **kw):
    """one change to the sliding seam's own contact (canal_seam_slide)."""
    def fn(m):
        assert m.set_contact('canal_seam_slide', **kw) == 1
        m.log.append(note)
    return fn


BUILDS = (('L93_pen5_paperP1_seamslide', 'L68_pen5_paperP1', seam_slide),
          ('L93_springs_newline_rhoi0_seamslide', 'L87_springs_newline_rhoi0', seam_slide),
          ('L93_newline_vwyeoh_rhoi0_seamslide', 'L85_newline_vwyeoh_rhoi0', seam_slide),
          # 13:45, the user: stop them and retry on P1. The Ogden runs stalled at t 0.23 / 0.25 (steps ~1e-6, the seam held);
          # with the faithful walls the strip nodes left the 2 mm search radius and the seam came apart (32 mm). One change
          # each from L93_pen5_paperP1_seamslide, and both:
          ('L94_pen5_paperP1_seamslide_rad10', 'L93_pen5_paperP1_seamslide',
           slide_set('the sliding seam: search_radius 2 -> 10 mm (nodes cannot escape it)', search_radius=10)),
          ('L94_pen5_paperP1_seamslide_pen50', 'L93_pen5_paperP1_seamslide',
           slide_set('the sliding seam: penalty 5 -> 50 (holds harder)', penalty=50)),
          ('L94_pen5_paperP1_seamslide_rad10_pen50', 'L93_pen5_paperP1_seamslide',
           slide_set('the sliding seam: search_radius 2 -> 10 mm and penalty 5 -> 50', search_radius=10, penalty=50)))
if __name__ == '__main__':
    if '--gaps' in sys.argv:
        m = Model7(os.path.join(RUNS, 'L87_springs_newline_rhoi0', 'L87_springs_newline_rhoi0.feb'))
        seam, strip, g = strip_gaps(m)
        print(f'{len(seam)} seam nodes, {len(strip)} strip facets, {len(g)} strip nodes; rest gap median {np.median(g):.3f}, '
              f'p10 {np.percentile(g, 10):.3f}, p90 {np.percentile(g, 90):.3f}, range {g.min():.3f}-{g.max():.3f} mm')
        sys.exit()
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        print(name, ':', m.log[-1][:160])
        emit(name, m)
