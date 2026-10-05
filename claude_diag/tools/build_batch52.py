"""Batch 52 (2026-09-26 ~09:10): the strongest connective-tissue-only change set on each line, adding a lower apex.

Springs line: P-arcus cut + side sphincters 50 % reached +44.2 deg at t = 1 (L40_springs_noparcus_sphside50); CL/USL at 50 %
added +6 deg at t 0.5 to P-arcus 10 % (L43_springs_parcus10_clusl50, running). Lofts line: the distal P-arcus cut + side
sphincters 50 % (L41_lofts_parcusdist0_sphside50) is the furthest lofts variant that does not crawl. One change each:
  L45_springs_noparcus_sphside50_clusl50   L40_springs_noparcus_sphside50 + every CL/USL connector set at 50 %
  L45_lofts_parcusdist0_sphside50_clusl50  L41_lofts_parcusdist0_sphside50 + the CL-L/R, USL-L/R lofts at 50 % (c1 and k
                                           x 0.5, as L30_lofts_clusl50_LA10kPa)
All NOT IN SOURCE (connective tissue).
usage: py -3.10 build_batch52.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch50 import clusl50

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
LOFTS = ('CL-L_fan', 'CL-R_fan', 'USL-L_fan', 'USL-R_fan')


def clusl_lofts50(m):
    doms = {d.get('name'): d.get('mat') for d in m.root.find('MeshDomains')}
    assert len({doms[l] for l in LOFTS}) == len(LOFTS), 'the CL/USL lofts share a material'
    mats = {x.get('name'): x for x in m.root.find('Material')}
    for loft in LOFTS:
        mat = mats[doms[loft]]
        assert mat.get('type') == 'Ogden', mat.get('type')
        assert all(float(mat.find(f'c{i}').text) == 0 for i in range(2, 7) if mat.find(f'c{i}') is not None)
        old = {k: mat.find(k).text for k in ('c1', 'k')}
        for k in ('c1', 'k'):
            mat.find(k).text = '%.6g' % (float(old[k]) * 0.5)
        m.log.append(f'NOT IN SOURCE (connective tissue): {loft} loft ({mat.get("name")}) at 50 %: c1 {old["c1"]} -> '
                     f'{mat.find("c1").text}, k {old["k"]} -> {mat.find("k").text} (m1 unchanged)')


BUILDS = (('L45_springs_noparcus_sphside50_clusl50', 'L40_springs_noparcus_sphside50', clusl50),
          ('L45_lofts_parcusdist0_sphside50_clusl50', 'L41_lofts_parcusdist0_sphside50', clusl_lofts50))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
