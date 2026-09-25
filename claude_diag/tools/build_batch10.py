"""Batch 10 (2026-09-23): one change each on top of TRUSS_fixed (= T3f_k_snn0)."""
import numpy as np

from variants2 import Model2, BEST, emit
from abq_surf import surface_facets, key, normal
from febmodel import Feb


def peb_bottom_facets():
    """Abaqus Surf-PVW-front-Peb-bottom facets absent from FEBio SlidingElastic1Primary, as FEBio node ids,
    ordered so the facet normal points out of the owning hex (FEBio convention)."""
    f = Feb(BEST)
    X = f.nodes
    have = {key(np.array([X[i] for i in c])) for _, c in f.surfaces['SlidingElastic1Primary']}
    ids = np.array(sorted(X))
    P = np.array([X[i] for i in ids])
    hexes = [c for et, d in f.elem_blocks.values() if et == 'hex8' for c in d.values()]
    owner = {}
    for c in hexes:
        for n in c:
            owner.setdefault(n, []).append(c)
    out = []
    for e, face, Q in surface_facets('Surf-PVW-front-Peb-bottom', 'VW-PeB'):
        if key(Q) in have:
            continue
        fid = []
        for q in Q:
            d = np.linalg.norm(P - q, axis=1)
            k = int(np.argmin(d))
            assert d[k] < 1e-6
            fid.append(int(ids[k]))
        own = [c for c in owner[fid[0]] if set(fid) <= set(c)]
        assert len(own) == 1, own
        cen = np.array([X[i] for i in own[0]]).mean(0)
        Fq = np.array([X[i] for i in fid])
        if normal(Fq) @ (Fq.mean(0) - cen) < 0:
            fid = fid[::-1]
        out.append(fid)
    return out


if __name__ == '__main__':
    m = Model2(BEST); m.add_pvw_la_contact(); emit('P1_pvwla', m)
    m = Model2(BEST); m.add_la_pressure(); emit('L1_laload', m)
    m = Model2(BEST); m.add_vw_mid_xsymm(); emit('B1_vwmid', m)
    m = Model2(BEST); m.set_connector_extend('constant'); emit('X1_extconst', m)
    fac = peb_bottom_facets(); print('PeB-bottom facets', len(fac))
    m = Model2(BEST); m.add_peb_bottom_to_avw_contact(fac); emit('C1_pebavw', m)
    m = Model2(BEST); m.set_diagnostic_output(); emit('G0_diag', m)
