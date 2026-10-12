"""Batch 243 (2026-10-09, the user's new goal: "We're no longer trying to match the numbers from the paper ... get about the same
range of motion ... make sure the movement of our model isn't too restrictive"; "Range is just that the Bp and Ba can get at
least to the points that they got to in the paper. If the paper has apical location check for that too. Healthy case should
stay unprolapsed."). The new standing settings: the remeshed PM, the cervix x3 of the wall, every support and wall at the
source value (no CL / Parcus fit scales); Load-LA x1/3, sphincter side x0.25, arcus chain mass x10, ledges x0.1 kept.
The paper's cases (Luo et al. 2015, 140 cm H2O; impairment = the healthy stiffness x (1 - impairment)):
  healthy  none                                                        Ba, Bp < 0
  P1       levator 20 %, apical 30 %, posterior 85 %                   Bp >= +4 (build_batch107 P1)
  P2       levator 60 %, apical 60 %, posterior 85 %                   Bp >= +9 (build_batch107 P2)
  F7A      levator 50 %, apical 80 %, anterior 50 %, posterior 90 %    rectocele, Bp >= +10 (Fig. 7A: Ba -1, Bp +1 cm)
  F8A      levator 50 %, anterior 90 %, posterior 50 %                 cystocele, Ba >= +10 (Fig. 8A: Ba +1, Bp -1 cm)
  F6a0     levator 20 %, posterior 85 %                                apex moves >= 0.18 mm / cm H2O x 140 = 25 mm
  F6a60    levator 20 %, apical 60 %, posterior 85 %                   apex moves >= 0.25 x 140 = 35 mm
           (Fig. 6, also anterior 0.17 / 0.21 -> 24 / 29 mm, posterior 0.29 / 0.36 -> 41 / 50 mm of motion)
Anterior support = the AVW-Para springs (the paper: "the origins of anterior vaginal support lay on the ATFP").
usage: py -3.10 build_batch243.py BATCH TAG [key=S ...] [--cases ...] [--base RUN]   (as build_batch145)"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_batch145 as b145  # noqa: E402
from build_batch44 import spring_scale  # noqa: E402
from build_batch50 import CLUSL  # noqa: E402
from build_batch107 import la_healthy_pct, LABEL  # noqa: E402


def case(la, apical=1.0, anterior=1.0, posterior=1.0):
    def fn(m):
        la_healthy_pct(m, la)
        if apical != 1:
            for s in CLUSL:
                spring_scale(m, s, apical, LABEL + ', apical')
        if anterior != 1:
            for s in ('AVW-Para-L_conn', 'AVW-Para-R_conn'):
                spring_scale(m, s, anterior, LABEL + ', anterior support')
        if posterior != 1:
            spring_scale(m, 'Parcus_conn', posterior, LABEL + ', posterior support')
    return fn


b145.CASES.update({'F7A': case(50, apical=0.2, anterior=0.5, posterior=0.1),
                   'F8A': case(50, anterior=0.1, posterior=0.5),
                   'F6a0': case(80, posterior=0.15),
                   'F6a60': case(80, apical=0.4, posterior=0.15),
                   # 2026-10-09 evening (the user: reach the ranges "maybe with different impairments"): the heaviest
                   # combination, levator 60 %, apical 80 %, anterior 50 %, posterior 95 % impaired
                   'MAX': case(40, apical=0.2, anterior=0.5, posterior=0.05)})

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('batch'); ap.add_argument('tag'); ap.add_argument('kv', nargs='*')
    ap.add_argument('--cases', default='healthy,P1,P2,F7A,F8A,F6a0,F6a60')
    ap.add_argument('--base', default='L228_ext5sl01_pmc_pmd_sph025')
    a = ap.parse_args()
    sc = {k: float(v) for k, v in (x.split('=') for x in a.kv)}
    print(' '.join(b145.build(a.batch, a.tag, sc, a.cases.split(','), a.base)))
