"""Batch 42 (2026-09-25, user-chosen): the PVW at 50 % stiffness on bases with a better-held apex, and the USL-L tie.

Batch 41 showed a softer PVW moves both walls forward (Vagina_PVW x 0.5 on L31: mean u_y PVW +13.6, AVW +11.3 mm at t 0.72,
against -2.4 / +0.4), but on L31's weak apex (CL/USL 30 %, sphincter connectors 50 %) the runs crawled: the apex snapped
(12 m/s, PVW x 0.5) or the USL-L loft flapped (PVW x 0.2). All at Load-LA 0.01, the other four 0.014; NOT IN SOURCE:
  L33_pvw50_clusl50         L30_lofts_clusl50_LA10kPa (CL/USL 50 %, sphincters 100 %) + Vagina_PVW x 0.5
  L33_pvw50_clusl100        L30_lofts_LA10kPa (CL/USL 100 %) + Vagina_PVW x 0.5
  L33_pvw50_clusl50_usltie  L33_pvw50_clusl50 + the USL-L loft on shared nodes like USL-R: its two tied-node-on-facet ties
                            (to _PickedSet64 = PVW, to _PickedSet346 = cervix) replaced by merging each tied tissue node with
                            the loft node it sits on (a loft convention; NOT IN SOURCE either way: the source has connectors)
usage: py -3.10 build_batch42.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

import numpy as np

from variants6 import Model6, JOBS, emit
from build_batch41 import scale_ogden  # noqa: E402  (build_batch41 builds only the names in sys.argv)
from febmodel import pt_tri_dist  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
USL_L_TIES = ('USL_L_fan__PickedSet64_tnof', 'USL_L_fan__PickedSet346_tnof')


def merge_tnof(m, contacts, loft='USL-L_fan', tol=1e-4):
    mesh = m.root.find('Mesh')
    X = m.nodes()
    surf = {s.get('name'): s for s in mesh.findall('Surface')}
    pairs = {p.get('name'): p for p in mesh.findall('SurfacePair')}
    blk = next(b for b in mesh.findall('Elements') if b.get('name') == loft)
    loft_nodes = {int(v) for e in blk for v in e.text.split(',')}
    remap, worst = {}, 0.0
    con = m.root.find('Contact')
    for cname in contacts:
        c = next(c for c in con if c.get('name') == cname)
        sp = pairs[c.get('surface_pair')]
        prim, sec = surf[sp.find('primary').text], surf[sp.find('secondary').text]
        pn = {int(v) for f in prim for v in f.text.split(',')}
        sn = {int(v) for f in sec for v in f.text.split(',')}
        assert sn <= loft_nodes, f'{cname}: secondary is not on {loft}'
        L = np.array(sorted(sn))
        LX = np.array([X[i] for i in L])
        maxd = float(c.find('max_distance').text)
        tris = [tuple(int(v) for v in f_.text.split(',')) for f_ in sec]
        T = [np.array([X[i] for i in t[:3]]) for t in tris]
        cen = np.array([t.mean(axis=0) for t in T])
        n_tied = 0
        for p in sorted(pn):
            x = np.array(X[p])
            d = np.linalg.norm(LX - x, axis=1)
            j = int(np.argmin(d))
            if d[j] >= tol:
                # the tie only takes primary nodes within max_distance of the loft facets: none may be left out
                if d[j] < maxd + 5.0:
                    ii = np.argsort(np.linalg.norm(cen - x, axis=1))[:20]
                    df = min(pt_tri_dist(x, *T[i]) for i in ii)
                    assert df > maxd, f'{cname}: tissue node {p} is {df:.3g} mm from the loft facets (tied, not merged)'
                continue
            worst = max(worst, float(d[j]))
            l = int(L[j])
            assert remap.get(l, p) == p, f'loft node {l} would map to two tissue nodes'
            remap[l] = p
            n_tied += 1
        assert n_tied > 0, cname
        con.remove(c)
        mesh.remove(sp)
        for s in (prim, sec):
            mesh.remove(s)
    # the re-pointed loft nodes must not be used anywhere else
    used = set()
    for tag in ('NodeSet',):
        for s in mesh.findall(tag):
            if s.text:
                used |= {int(v) for v in s.text.replace('\n', ',').split(',') if v.strip()}
    for s in mesh.findall('Surface'):
        used |= {int(v) for f in s for v in f.text.split(',')}
    for ds in mesh.findall('DiscreteSet'):
        used |= {int(v) for e in ds.findall('delem') for v in e.text.split(',')}
    clash = sorted(set(remap) & used)
    assert not clash, f'loft nodes used elsewhere: {clash[:10]}'
    for e in blk:
        e.text = ','.join(str(remap.get(int(v), int(v))) for v in e.text.split(','))
    m.log.append(f'{loft} on shared nodes (like USL-R): the tied-node-on-facet ties {list(contacts)} removed (with their '
                 f'surface pairs and surfaces) and the {len(remap)} loft nodes they tied re-pointed to the tissue nodes they '
                 f'sat on (max distance {worst:.2g} mm). A loft convention, NOT IN SOURCE (the source has connectors)')


BUILDS = (('L33_pvw50_clusl50', 'L30_lofts_clusl50_LA10kPa', False),
          ('L33_pvw50_clusl100', 'L30_lofts_LA10kPa', False),
          ('L33_pvw50_clusl50_usltie', 'L33_pvw50_clusl50', True))
for name, base, tie in BUILDS:
    if WANT and name not in WANT:
        continue
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(os.path.join(RUNS, base, base + '.feb'))
    if tie:
        merge_tnof(m, USL_L_TIES)
    else:
        scale_ogden(m, 'Vagina_PVW', 0.5)
    emit(name, m)
