"""Geometry facts for the Abaqus Int-PVW-LA contact (Surf-VW-PVW-back vs _PickedSurf479 = LA SPOS).

1. Which FEBio <Surface> is the Abaqus PVW back face (coordinate match), and is its normal outward?
2. Is the FEBio LA connectivity oriented like the Abaqus LA (so FEBio top face normal = Abaqus SPOS)?
3. Which side of the LA the PVW back face is on, and how far (signed, along the LA normal).
usage: py -3.10 pvw_la_geom.py model.feb
"""
import sys

import numpy as np

from abq_surf import surface_facets, part_elements, key, normal
from febmodel import Feb, pt_tri_dist
from inpmap import part_nodes

feb = Feb(sys.argv[1])
X = feb.nodes


def feb_facets(name):
    return [np.array([X[i] for i in c]) for _, c in feb.surfaces[name]]


# ---- 1. PVW back face -------------------------------------------------------
abq = {s: surface_facets(s, 'VW-PeB') for s in
       ('Surf-VW-PVW-back', 'Surf-PVW-front-Peb-bottom', 'Surf-VW-AVW-back', 'Surf-VW-AVW-front')}
abq_keys = {s: {key(P): normal(P) for _, _, P in f} for s, f in abq.items()}
print('Abaqus surfaces (VW-PeB):', {s: len(v) for s, v in abq_keys.items()})
for fs in feb.surfaces:
    if fs.startswith(('Load-', 'SlidingElastic')):
        F = feb_facets(fs)
        fk = {key(P): normal(P) for P in F}
        for s, ak in abq_keys.items():
            common = set(fk) & set(ak)
            if common:
                same = sum(float(fk[k] @ ak[k]) > 0 for k in common)
                print(f'  FEBio {fs:26s} ({len(fk)}) ~ Abaqus {s:26s} ({len(ak)}): '
                      f'{len(common)} common facets, normal same/flipped {same}/{len(common) - same}')

# ---- 2. LA orientation ------------------------------------------------------
LX = part_nodes('LA-new')
LE = part_elements('LA-new')
abq_la = {key([LX[i] for i in c]): normal(np.array([LX[i] for i in c])) for c in LE.values()}
feb_la = {}
for name in ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri'):
    for c in feb.elem_blocks[name][1].values():
        P = np.array([X[i] for i in c])
        feb_la[key(P)] = normal(P)
common = set(abq_la) & set(feb_la)
same = sum(float(abq_la[k] @ feb_la[k]) > 0 for k in common)
print(f'LA elements: Abaqus {len(abq_la)}, FEBio {len(feb_la)}, matched {len(common)}, '
      f'normal same/flipped {same}/{len(common) - same}')

# ---- 3. PVW back face relative to the LA mid-surface -------------------------
tris, tn = [], []
for c in LE.values():
    P = np.array([LX[i] for i in c])
    n = normal(P)
    for t in ([0, 1, 2], [0, 2, 3]) if len(c) == 4 else ([0, 1, 2],):
        tris.append(P[t])
        tn.append(n)
tris = np.array(tris)
tn = np.array(tn)
cen = tris.mean(axis=1)
pv_nodes = {tuple(np.round(p, 6)): p for _, _, P in abq['Surf-VW-PVW-back'] for p in P}
res = []
for p in pv_nodes.values():
    near = np.argsort(np.linalg.norm(cen - p, axis=1))[:40]
    d = [pt_tri_dist(p, *tris[i]) for i in near]
    k = near[int(np.argmin(d))]
    s = float((p - cen[k]) @ tn[k])
    res.append((min(d), s))
res = np.array(res)
close = res[res[:, 0] < 6]
print(f'PVW back nodes: {len(res)}; distance to LA mid-surface min {res[:, 0].min():.3f} mm; '
      f'nodes within 6 mm: {len(close)}, of which on +n (SPOS) side: {(close[:, 1] > 0).sum()}')
for lim in (2.5, 3, 3.5, 4, 5, 6):
    print(f'   within {lim} mm of LA mid-surface: {(res[:, 0] < lim).sum()}')
