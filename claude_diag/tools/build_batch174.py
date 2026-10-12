"""Batch 174 (2026-10-05, the user: "in all future models the LA_sphincter_side_conn_mat need to be less stiff. Start with 50%
less stiff ... it should allow Bp to increase" -> "run both").
LA_sphincter_side_conn: 16 springs between the levator (LA_PCMPRM) and the perineal body (_PickedSet66), source law
LA_sphincter_side_conn_mat; its <scale> 1 -> 0.5 (NOT IN SOURCE), one change on two P1 bases where Bp is short:
  L174_selfall_walls035_sphside05_P1   L169_wrapfull_endspr_selfall_walls035_P1 (tube, walls x0.35, damped, end springs,
                                       whole-boundary self-contact) + the change
  L174_avwcx5_pmdamp_sphside05_P1      L170_avwcx5_pvwband3_pmdamp_P1 (non-tube, contact-only seam, AVW + cervix c1 x5, PVW
                                       3 mm band, wall + PM damping) + the change
usage: py -3.10 build_batch174.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
MAT = 'LA_sphincter_side_conn_mat'
FACTOR = 0.5
BUILDS = (('L174_selfall_walls035_sphside05_P1', 'L169_wrapfull_endspr_selfall_walls035_P1'),
          ('L174_avwcx5_pmdamp_sphside05_P1', 'L170_avwcx5_pvwband3_pmdamp_P1'))


def soften(t):
    i = t.index(f'name="{MAT}" type="nonlinear spring">')
    j = t.index('</discrete_material>', i)
    blk, k = re.subn(r'<scale>([^<]+)</scale>', lambda m: f'<scale>{float(m.group(1)) * FACTOR:g}</scale>', t[i:j])
    assert k == 1
    old = re.search(r'<scale>([^<]+)</scale>', t[i:j]).group(1)
    return t[:i] + blk + t[j:], (f'NOT IN SOURCE (the user\'s rule from 2026-10-05): {MAT} scale {old} -> '
                                 f'{float(old) * FACTOR:g} (the 16 LA_sphincter_side_conn springs, levator <-> perineal body)')


if __name__ == '__main__':
    for name, src in BUILDS:
        p = os.path.join(RUNS, src, src + '.feb')
        t, msg = soften(open(p, encoding='utf-8').read())
        write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch174.py)\n{msg}\n')
