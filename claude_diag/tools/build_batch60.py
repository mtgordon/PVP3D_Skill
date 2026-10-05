"""Batch 60 (2026-09-26 ~10:45): the paper's Figure 3 case on both lines.

The user's target is Figure 3C of Luo et al. 2015 (J Biomech, "A multi-compartment 3-D finite element model of rectocele
and its interaction with cystocele"; the PDF is in the user's Downloads; text in scratch_2026-09-26/luo_rectocele_paper.txt):
"50% levator impairment, 70% apical impairment and 85% posterior support impairment, and with no anterior support
impairment", at 140 cm H2O (= 0.0137 MPa, the source's 0.014). Impairment = reduced stiffness; 90 % = "totally detached".
In our model: posterior support (PPS) = the P-arcus connectors (LA-Y-parcus-*); apical = CL/USL; the levator is already
LA_Yamada50% in the source; the sphincter (PeB-LA) connectors are not impaired in the paper. One change each, NOT IN SOURCE
(connective tissue), Load-LA kept at 1/3 (the user's setting):
  L53_springs_pc10_cu30               L40_springs_parcus10 (P-arcus 10 %) + every CL/USL connector set at 30 %
  L53_lofts_pc10_cu30_avwparaconn     L46_lofts_parcus10_avwparaconn + the CL/USL lofts at 30 % (c1, k x 0.3)
usage: py -3.10 build_batch60.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import spring_scale
from build_batch50 import CLUSL

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CT = 'NOT IN SOURCE (connective tissue)'
LOFTS = ('CL-L_fan', 'CL-R_fan', 'USL-L_fan', 'USL-R_fan')


def clusl30(m):
    for s in CLUSL:
        spring_scale(m, s, 0.3, CT)


def clusl_lofts30(m):
    doms = {d.get('name'): d.get('mat') for d in m.root.find('MeshDomains')}
    assert len({doms[l] for l in LOFTS}) == len(LOFTS)
    mats = {x.get('name'): x for x in m.root.find('Material')}
    for loft in LOFTS:
        mat = mats[doms[loft]]
        old = {k: mat.find(k).text for k in ('c1', 'k')}
        for k in ('c1', 'k'):
            mat.find(k).text = '%.6g' % (float(old[k]) * 0.3)
        m.log.append(f'{CT}: {loft} loft ({mat.get("name")}) at 30 %: c1 {old["c1"]} -> {mat.find("c1").text}, '
                     f'k {old["k"]} -> {mat.find("k").text} (m1 unchanged)')


BUILDS = (('L53_springs_pc10_cu30', 'L40_springs_parcus10', clusl30),
          ('L53_lofts_pc10_cu30_avwparaconn', 'L46_lofts_parcus10_avwparaconn', clusl_lofts30))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
