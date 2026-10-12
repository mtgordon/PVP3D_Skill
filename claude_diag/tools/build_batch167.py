"""Batch 167 (2026-10-05, the user: "go back to the non-tube version and try P1 and P2 with the stiff band to see if that helps
produce the paper results"; base L149 with the seam springs out, the wall damping, P1 + P2 only).
The stiff band (build_batch100.band, 2026-09-28: L102_newline_vwyeoh_rhoi0_band3x3 on L85, healthy only) was the way the
faithful (Yeoh) walls got through with a contact-only canal seam, before the seam springs (L91) and the tube: Yeoh c1 x3 in
the wall elements whose centroid lies within 3 mm of a lateral canal-seam node (~10 %), the rest faithful (NOT IN SOURCE).
From L149_seamspr_pm_outer (the newest non-tube base: L91 + seg_up 2 + the PM structure + PM_conn on the outer arc + walls_LA /
walls_PM):
  1. the canal seam springs removed (constraint canal_seam_springs; the seam is contact-only again, as in the source);
  2. the 3 mm band, c1 x3 (build_batch100.band(3, 3));
  3. the wall damping (build_batch165.wall_damping on the AVW, cervix and PVW: C 200/s from t = 0; NOT IN SOURCE);
  4. the paper's P1 / P2 (build_batch145).
  L167_pm_outer_band3x3_damp              steps 1-3 (healthy; not run)
  L167_pm_outer_band3x3_damp_{P1,P2}      + the cases
usage: py -3.10 build_batch167.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch100 import band  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch165 import wall_damping  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L149_seamspr_pm_outer'
TMP = 'L167_pm_outer_band3x3_tmp'
NAME = 'L167_pm_outer_band3x3_damp'
WALLS = ('_PickedSet347', '_PickedSet346', '_PickedSet64')

if __name__ == '__main__':
    for nm in (TMP, NAME):
        assert not os.path.exists(os.path.join(RUNS, nm)), nm
    m = Model7(os.path.join(RUNS, SRC, SRC + '.feb'))
    cons = m.root.find('Constraints')
    seam = [c for c in cons if c.get('name') == 'canal_seam_springs']
    assert len(seam) == 1
    n = len(seam[0].findall('linear_constraint'))
    cons.remove(seam[0])
    m.log.append(f'as the source: the canal seam springs removed (constraint canal_seam_springs, {n} linear constraints); the '
                 f'seam is contact-only (SlidingElastic1)')
    band(3, 3)(m)
    emit(TMP, m)
    p = os.path.join(RUNS, TMP, TMP + '.feb')
    t, msg = wall_damping(open(p, encoding='utf-8').read(), WALLS)
    notes = open(p + '.changes.txt', encoding='utf-8').read()
    write(NAME, t, f'{notes}{msg}\n(tools/build_batch167.py; base {SRC})\n')
    import shutil
    shutil.rmtree(os.path.join(RUNS, TMP))
    print(' '.join(build(167, NAME[len('L167_'):], {}, ['P1', 'P2'], NAME)))
