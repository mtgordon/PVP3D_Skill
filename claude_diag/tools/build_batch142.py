"""Batch 142 (2026-10-01, the user: "Can you wrap the edges all the way to the end of the vaginal wall?" -> option 2: wrap to
the end at the perineal-body end and let the curve flare with the PVW; leave the top as is; "try and run it. Fix it if you
find issues and try running again").
As L141_thin275_narrow3_wrapcontact (L137_tube_thin275_narrow3_flare_pmin + the curves' inner faces in the lumen contact), with:
  * the curves carried to the last seam station at the perineal-body end on each side (build_batch135.EXTEND_DISTAL: those
    stations get PVW partners further along the PVW lumen edge, PeB-shared nodes allowed; the top end is unchanged: the
    cervix already joins the PVW there);
  * the flare on the whole wall end: the AVW, cervix and PVW all widen back to full width over the last 3.4 mm before the
    PVW / perineal-body junction (wall_reshape flare_doms = all three), so the curves, built on the walls' columns, flare
    with them (before: the PVW only, and the curves stopped where the flare began);
  * the PM_conn family on the PM's inner arc now includes the middle spring (PM-LA-x%stiff_mat): 27 springs, one mean scale.
  L142_thin275_narrow3_wrapend
usage: py -3.10 build_batch142.py"""
import os
import re
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch90 import seam_nodes  # noqa: E402
from build_batch118 import fold  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM  # noqa: E402
from build_batch132 import seg_up2  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402
import wall_reshape  # noqa: E402
import build_batch135 as b135  # noqa: E402
import build_batch136 as b136  # noqa: E402
import build_batch141 as b141  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
NAME = os.environ.get('B142_NAME', 'L142_thin275_narrow3_wrapend')
PRE, MID = NAME + '_tmppre', NAME + '_tmpmid'
TOOLS = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    os.makedirs(SCRATCH, exist_ok=True)
    for nm in (NAME, PRE, MID):
        assert not os.path.exists(os.path.join(RUNS, nm)), f'{nm} exists; not overwriting'
    orig = os.path.join(RUNS, BASE, BASE + '.feb')
    wall_reshape.NARROW = 3.0
    spec = b136.flare_spec(orig)
    m0 = Model7(orig)
    rep = wall_reshape.reshape(m0, flare=spec, flare_doms=(tube.AVW, tube.CX, tube.PVW))
    mid = os.path.join(SCRATCH, NAME + '_reshaped.feb')
    m0.write(mid)
    m = Model113(mid)
    m.log.append(f'{b135.LABEL}: wall_reshape (narrowed 3 mm; the whole wall end flared back to full width over '
                 f'{[round(v[1], 1) for v in spec.values()]} mm before the PeB junction; the PeB not narrowed): {rep}')
    b135.NAME = NAME
    b135.EXTEND_DISTAL, b135.EXTENDED = True, {}
    jh = b135.wrap2(m, orig)
    m.log.append(f'the curves carried to the end at the PeB end: extra station partners (AVW seam node -> PVW lumen-edge '
                 f'node) {b135.EXTENDED}')
    print(m.log[-2][:400]); print(m.log[-1])
    seg_up2(m)
    step2(INNER + OUTER_BOTTOM)(m)
    inside = b135.walls_la(m)
    print(m.log[-1][:400])
    emit(PRE, m)
    b136.run_main(os.path.join(TOOLS, 'build_batch129.py'), {'BASE': PRE, 'NAME': MID})
    b136.run_main(os.path.join(TOOLS, 'build_batch129b.py'), {'SRC': MID, 'OLD': PRE, 'NAME': NAME + '_tmpavg'})
    # the curves' inner faces into the lumen contact (batch 141)
    seam, _, _ = seam_nodes(Model7(orig))
    rp, _ = fold(Model7(orig))
    inv = {a: p for p, a in rp.items()}
    inv.update(b135.EXTENDED)
    src = os.path.join(RUNS, NAME + '_tmpavg', NAME + '_tmpavg.feb')
    f, t = Feb(src), open(src, encoding='utf-8').read()
    C = np.array([0.5 * (np.array(f.nodes[a]) + np.array(f.nodes[inv[a]])) for a in seam if a in inv])
    fp = b141.inner_faces(f, 'canal_wrap_hex_pvw', C)
    fa = b141.inner_faces(f, 'canal_wrap_hex_avw', C)
    t2 = b141.add_faces(b141.add_faces(t, 'SlidingElastic1Primary', fp), 'SlidingElastic1Secondary', fa)
    out = os.path.join(RUNS, NAME)
    os.makedirs(out)
    open(os.path.join(out, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
    notes = [open(os.path.join(RUNS, nm, nm + '.feb.changes.txt'), encoding='utf-8').read() for nm in (PRE, MID, NAME + '_tmpavg')]
    with open(os.path.join(out, NAME + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
        fo.write(f'base: {orig} (tools/build_batch142.py)\n' + '\n'.join(notes)
                 + f'\nthe curves\' inner faces in the lumen contact: {len(fp)} PVW-half faces -> SlidingElastic1Primary, '
                   f'{len(fa)} AVW-half faces -> SlidingElastic1Secondary\n')
    for nm in (PRE, MID, NAME + '_tmpavg'):
        shutil.move(os.path.join(RUNS, nm), os.path.join(SCRATCH, nm))
    print('inner faces', len(fp), len(fa), '; wall nodes inside the LA at rest', len(inside))
    print('wrote', os.path.join(out, NAME + '.feb'))
