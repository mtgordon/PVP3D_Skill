"""Batch 168 (2026-10-05, the user: "do you think that you could get the stiffer seams (maybe expanding the stiffness further
in) to get close to the P1, P2 numbers?" -> "build the three band P1 runs").
On the band base L167_pm_outer_band3x3_damp (L149 with the seam springs out (contact-only seam, as the source), Yeoh c1 x3
within 3 mm of the seam, the wall damping): stiffening moves Ba (the AVW and cervix descend; the PVW barely moves), so:
  1. L168_avwcx3_pvwband3_P1    AVW + cervix c1 x3 everywhere (their c1 maps set uniform), the PVW keeps the 3 mm band
  2. L168_avwcx5_pvwband3_P1    AVW + cervix c1 x5 everywhere (the top of the user's bounds), the PVW keeps the 3 mm band
  3. L168_band6x3_P1            the band widened to 6 mm (c1 x3) in all three walls (build_batch100.band(6, 3) on L149
                                with the seam springs out, as build_batch167), the wall damping
each + the paper's P1 (build_batch145). All NOT FAITHFUL / NOT IN SOURCE (c1 scaled; c2 and k kept, as the band).
usage: py -3.10 build_batch168.py"""
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch100 import band  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch165 import wall_damping  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BAND = 'L167_pm_outer_band3x3_damp'
SRC = 'L149_seamspr_pm_outer'
C1 = 0.005964127
WALLS = ('_PickedSet347', '_PickedSet346', '_PickedSet64')


def uniform(t, maps, factor):
    for mp in maps:
        i = t.index(f'<ElementData name="{mp}"')
        j = t.index('</ElementData>', i)
        blk, k = re.subn(r'(<elem lid="\d+">)[^<]+(</elem>)', lambda m: f'{m.group(1)}{C1 * factor:.7g}{m.group(2)}', t[i:j])
        assert k > 100, (mp, k)
        t = t[:i] + blk + t[j:]
    return t


if __name__ == '__main__':
    p = os.path.join(RUNS, BAND, BAND + '.feb')
    t0 = open(p, encoding='utf-8').read()
    notes = open(p + '.changes.txt', encoding='utf-8').read()
    for f in (3, 5):
        name = f'L168_avwcx{f}_pvwband3'
        t = uniform(t0, ('Vagina_AVW_c1map', 'Vagina_Cervix_c1map'), f)
        write(name, t, f'base: claude_diag/runs/{BAND}/{BAND}.feb (tools/build_batch168.py)\n{notes}NOT FAITHFUL: Vagina_AVW and '
                       f'Vagina_Cervix Yeoh c1 x{f} ({C1 * f:.7g}) in every element (their c1 maps set uniform; c2, k kept); '
                       f'Vagina_PVW keeps the 3 mm band\n')
        print(' '.join(build(168, name[len('L168_'):], {}, ['P1'], name)))
    # 3. the 6 mm band in all walls
    tmp, name = 'L168_band6x3_tmp', 'L168_band6x3'
    m = Model7(os.path.join(RUNS, SRC, SRC + '.feb'))
    cons = m.root.find('Constraints')
    seam = [c for c in cons if c.get('name') == 'canal_seam_springs']
    cons.remove(seam[0])
    m.log.append('as the source: the canal seam springs removed (constraint canal_seam_springs); the seam is contact-only')
    band(6, 3)(m)
    emit(tmp, m)
    q = os.path.join(RUNS, tmp, tmp + '.feb')
    t, msg = wall_damping(open(q, encoding='utf-8').read(), WALLS)
    write(name, t, open(q + '.changes.txt', encoding='utf-8').read() + msg + f'\n(tools/build_batch168.py; base {SRC})\n')
    shutil.rmtree(os.path.join(RUNS, tmp))
    print(' '.join(build(168, name[len('L168_'):], {}, ['P1'], name)))
