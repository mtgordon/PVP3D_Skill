"""Batch 75 (2026-09-27): two published cases of the source model, to compare FEBio with the Abaqus results in numbers.

Luo et al. 2015 (J Biomech, the source model's paper; text in scratch_2026-09-26/luo_rectocele_paper.txt), levator series
at 140 cm H2O (= the source's 0.014 MPa) with 85 % posterior support impairment and normal anterior support: "the Bp values
change from 0.4 cm to 0.9 cm when levator and apical impairments were increased from 20% and 30% to 60% and 60%", and all
Ba values are above the hymen (< 0). Bp / Ba = the posterior / anterior wall points 3 cm above the hymenal ring
(tools/pop_q.py measures them against PM_Plane, the source's hymenal plane).
Impairment = reduced tensile stiffness (the paper). Mapping (as on 2026-09-26, approximate: the paper's 0 % baselines were
calibrated and not published; the source .inp values are taken as the 0 % baseline):
  levator: the source LA material LA_Yamada50% is 50 % impaired; LA_YamadaX% = X/50 x LA_Yamada50% (the .inp's 10 / 50 /
           80 % tables are exact multiples), so the fitted Yeoh c1, c2 and k scale by X/50;
  apical: every CL/USL connector set; posterior support: the 26 P-arcus connectors (Parcus_conn).
A published case sets several parameters by definition (not a one-change variant):
  P1: levator 20 % (LA x 1.6), apical 30 % (CL/USL x 0.7), posterior 85 % (P-arcus x 0.15): paper Bp = +4 mm
  P2: levator 60 % (LA x 0.8), apical 60 % (CL/USL x 0.4), posterior 85 % (P-arcus x 0.15): paper Bp = +9 mm
On the fast base (L19_pen5, chain mass x10) and like for like (L21D_allconn_pen, the source chain mass):
  L68_pen5_paperP1, L68_pen5_paperP2, L68_d_paperP1, L68_d_paperP2
usage: py -3.10 build_batch75.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch44 import spring_scale
from build_batch50 import CLUSL

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
LABEL = 'the paper\'s case (Luo et al. 2015)'


def la_scale(m, s):
    mat = m._material('LA_Yamada50pct_Yeoh')
    old = {k: mat.find(k).text for k in ('c1', 'c2', 'k')}
    for k in ('c1', 'c2', 'k'):
        mat.find(k).text = '%.6g' % (float(old[k]) * s)
    m.log.append(f'{LABEL}: levator at {50 * s:g} % of the healthy Yamada stiffness (LA_Yamada{50 * s:g}%): '
                 f'LA_Yamada50pct_Yeoh c1, c2, k x {s:g}: {old} -> '
                 f'{ {k: mat.find(k).text for k in ("c1", "c2", "k")} }')


def case(la, clusl, parcus):
    def fn(m):
        la_scale(m, la)
        for s in CLUSL:
            spring_scale(m, s, clusl, LABEL + ', apical')
        spring_scale(m, 'Parcus_conn', parcus, LABEL + ', posterior support')
    return fn


P1 = case(1.6, 0.7, 0.15)
P2 = case(0.8, 0.4, 0.15)
BUILDS = (('L68_pen5_paperP1', 'L19_pen5', P1), ('L68_pen5_paperP2', 'L19_pen5', P2),
          ('L68_d_paperP1', 'L21D_allconn_pen', P1), ('L68_d_paperP2', 'L21D_allconn_pen', P2))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
