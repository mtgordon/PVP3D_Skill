"""Batch 107 (2026-09-28 ~18:55, the user's choice: "Paper cases with softer walls"): the paper's cases P1 / P2 (Luo et al.
2015; build_batch75.py's mapping) on the current springs line, with the line's Ogden walls and with the candidate softer
walls (Yeoh c1 x1.5 + dtol 0.01), to see how the wall law changes the trend across load cases before adopting it.

The current line's LA is the healthy law (the source's PCM-LA_Yamada100% table as Yeoh: c1 0.02347128, c2 0.03283498,
k 2; the material keeps the name LA_Yamada50pct_Yeoh), with Load-LA x 1/3 (the user's choice for the lines; kept here).
As in build_batch75.py the source's connector values are the 0 % baseline and the Yamada tables are exact multiples:
  P1: levator 20 % impaired = 80 % of healthy (LA c1, c2, k x 0.8), apical 30 % (every CL / USL spring set x 0.7),
      posterior 85 % (Parcus_conn x 0.15); the paper: Bp +4 mm below the hymen, Ba above it
  P2: levator 60 % = 40 % of healthy (x 0.4), apical 60 % (x 0.4), posterior 85 % (x 0.15); the paper: Bp +9 mm, Ba above
A published case sets several parameters by definition (not a one-change variant); each pair differs only in the walls.
  L107_springs_newline_rhoi0_paperP1 / _paperP2                 L87_springs_newline_rhoi0 (the springs line) + P1 / P2
  L107_newline_vwyeoh_rhoi0_c1x15_dtol01_paperP1 / _paperP2     L103_newline_vwyeoh_rhoi0_c1x15_dtol01 + P1 / P2
usage: py -3.10 build_batch107.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch44 import spring_scale
from build_batch50 import CLUSL

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
LABEL = 'the paper\'s case (Luo et al. 2015)'
LA100 = {'c1': 0.02347128, 'c2': 0.03283498, 'k': 2.0}


def la_healthy_pct(m, pct):
    mat = m._material('LA_Yamada50pct_Yeoh')
    old = {k: mat.find(k).text for k in ('c1', 'c2', 'k')}
    assert all(abs(float(old[k]) - v) < 1e-9 for k, v in LA100.items()), old
    for k in ('c1', 'c2', 'k'):
        mat.find(k).text = '%.7g' % (LA100[k] * pct / 100)
    m.log.append(f'{LABEL}: levator at {pct:g} % of the healthy Yamada stiffness: the healthy law (PCM-LA_Yamada100%, '
                 f'material LA_Yamada50pct_Yeoh) c1, c2, k x {pct / 100:g}: {old} -> '
                 f'{ {k: mat.find(k).text for k in ("c1", "c2", "k")} }')


def case(la_pct, clusl, parcus):
    def fn(m):
        la_healthy_pct(m, la_pct)
        for s in CLUSL:
            spring_scale(m, s, clusl, LABEL + ', apical')
        spring_scale(m, 'Parcus_conn', parcus, LABEL + ', posterior support')
    return fn


P1 = case(80, 0.7, 0.15)
P2 = case(40, 0.4, 0.15)
SOFT = 'L103_newline_vwyeoh_rhoi0_c1x15_dtol01'
BUILDS = (('L107_newline_vwyeoh_rhoi0_c1x15_dtol01_paperP1', SOFT, P1),
          ('L107_springs_newline_rhoi0_paperP1', 'L87_springs_newline_rhoi0', P1),
          ('L107_newline_vwyeoh_rhoi0_c1x15_dtol01_paperP2', SOFT, P2),
          ('L107_springs_newline_rhoi0_paperP2', 'L87_springs_newline_rhoi0', P2))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
