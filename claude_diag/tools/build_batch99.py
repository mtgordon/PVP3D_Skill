"""Batch 99 (2026-09-28 ~15:20, the user away ~4-6 h): how close to the source's wall softness the vaginal walls can get
before FEBio's implicit solver stalls (the user, 2026-09-28: "keep investigating how to get softer vaginal wall materials
without it causing problems").

The faithful (Yeoh) walls (c1 0.005964127, c2 0.01883547, k 1; AVW, PVW, cervix) crawl at the canal's contact-only side
seam from t ~0.78 on the springs line (L85_newline_vwyeoh_rhoi0). The converter's Ogden walls converge; they are ~3x
stiffer at 20 % uniaxial strain (initial shear modulus 0.060 vs 0.012 MPa) but 2x softer in equibiaxial stretch 1.5.
Raise only the small-strain term c1, keeping c2 and k: NOT FAITHFUL, a deliberate compromise ("the most faithful walls
that still converge"). Uniaxial nominal stress relative to the faithful law at stretch 1.1 / 1.2 / 1.5 (tools/wall_laws.py):
c1 x1.5 1.42 / 1.30 / 1.11, x2 1.85 / 1.60 / 1.21, x3 2.70 / 2.19 / 1.43 (the Ogden walls 4.25 / 2.99 / 1.11).
One change each from L85_newline_vwyeoh_rhoi0 (= the springs line L87_springs_newline_rhoi0 + the Yeoh walls):
  L99_newline_vwyeoh_rhoi0_c1x15    Vagina_AVW, Vagina_PVW, Vagina_Cervix: Yeoh c1 x1.5
  L99_newline_vwyeoh_rhoi0_c1x2     the three walls: c1 x2
  L99_newline_vwyeoh_rhoi0_c1x3     the three walls: c1 x3 (the bracket toward the Ogden walls)
  L99_newline_vwyeoh_rhoi0_avwc1x2  Vagina_AVW only: c1 x2 (the seam mostly opens: the AVW's edges lift off the PVW; the
                                    PVW and cervix stay faithful)
Built on demand to narrow the bracket (named only): L101_newline_vwyeoh_rhoi0_c1x25 (x2.5), L101_newline_vwyeoh_rhoi0_c1x125
(x1.25), L101_newline_vwyeoh_rhoi0_avwc1x3 (the AVW alone x3).
usage: py -3.10 build_batch99.py [NAME ...]   (default: batch 99; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
WALLS = ('Vagina_AVW', 'Vagina_PVW', 'Vagina_Cervix')
C1 = 0.005964127                      # the faithful Yeoh c1 (build_batch77.VW)


def c1x(f, walls=WALLS):
    def fn(m):
        for n in walls:
            mat = m._material(n)
            assert mat.get('type') == 'Yeoh' and abs(float(mat.find('c1').text) - C1) < 1e-9, n
            mat.find('c1').text = '%.7g' % (C1 * f)
        m.log.append(f'NOT FAITHFUL (a deliberate compromise): {", ".join(walls)} Yeoh c1 {C1} -> {C1 * f:.7g} (x{f}; '
                     f'c2 0.01883547 and k 1 kept): the small-strain stiffness raised toward the Ogden walls to test '
                     f'how soft the walls can be before the implicit solver stalls at the canal seam')
    return fn


BUILDS = (('L99_newline_vwyeoh_rhoi0_c1x15', c1x(1.5)),
          ('L99_newline_vwyeoh_rhoi0_c1x2', c1x(2)),
          ('L99_newline_vwyeoh_rhoi0_c1x3', c1x(3)),
          ('L99_newline_vwyeoh_rhoi0_avwc1x2', c1x(2, ('Vagina_AVW',))),
          # built on demand to narrow the bracket once batch 99 shows where the stall starts (same base, one change each)
          ('L101_newline_vwyeoh_rhoi0_c1x25', c1x(2.5)),
          ('L101_newline_vwyeoh_rhoi0_c1x125', c1x(1.25)),
          ('L101_newline_vwyeoh_rhoi0_avwc1x3', c1x(3, ('Vagina_AVW',))))
LATER = {'L101_newline_vwyeoh_rhoi0_c1x25', 'L101_newline_vwyeoh_rhoi0_c1x125', 'L101_newline_vwyeoh_rhoi0_avwc1x3'}
if __name__ == '__main__':
    base = 'L85_newline_vwyeoh_rhoi0'
    for name, fn in BUILDS:
        if (WANT and name not in WANT) or (not WANT and name in LATER):
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
