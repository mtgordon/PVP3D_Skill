"""Batch 136 (2026-10-01, the user: "build another model that flares out the pvw from where the curved edges stop to where
it joins with the pe_b instead of narrowing the pe_b"; and "Make all future models attach the PM_Plane": the PM_conn
springs on the PM's inner arc).
As L135_tube_thin275_narrow6 (build_batch135: 2.75 mm walls with their outer surfaces kept, the whole wall 6 mm narrower
on each side, the curved edges 2.75 mm thick in AVW / PVW halves, PVW_LA + walls_LA contacts, seg_up 2, the PM as a
structure), with two changes:
  * the PVW flares: its narrowing fades from full where the curved edges stop (the distal-most wrapped station on each
    side) to none at the PVW / perineal-body junction (smoothstep in the distance to the junction nodes), and the
    perineal body follows only the PVW's thinning at the shared face, not the narrowing (wall_reshape.reshape(flare=...));
  * the 26 PM_conn springs attached to the PM's inner arc (build_batch129), PM_conn_mat's elongation scaled by their mean
    length change (build_batch129b), as L129_tube_r30s_pm_conninner_avg.
  L136_tube_thin275_flare_pmin   built only (the user looks at it first)
usage: py -3.10 build_batch136.py"""
import os
import shutil
import sys
import xml.etree.ElementTree as ET  # noqa: F401

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch118 import fold  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM  # noqa: E402
from build_batch132 import seg_up2  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402
import wall_reshape  # noqa: E402
import build_batch135 as b135  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
NAME = 'L136_tube_thin275_flare_pmin'
PRE, MID = 'L136tmp_pre', 'L136tmp_mid'


def run_main(module_file, overrides):
    """run a builder script's __main__ block with its constants overridden"""
    import importlib
    mod = importlib.import_module(os.path.splitext(os.path.basename(module_file))[0])
    for k, v in overrides.items():
        setattr(mod, k, v)
    src = open(module_file, encoding='utf-8').read()
    main = src[src.rindex("if __name__ == '__main__':"):].replace("if __name__ == '__main__':", 'if True:', 1)
    exec(compile(main, module_file, 'exec'), mod.__dict__)


def flare_spec(orig):
    f = Feb(orig)
    X = {n: np.array(p, float) for n, p in f.nodes.items()}
    junction = sorted(set(f.domain_nodes(tube.PVW)) & set(f.domain_nodes(tube.PEB)))
    rep, _ = fold(Model7(orig))                      # PVW lumen-edge nodes that the curved edges attach to
    x0 = float(np.mean([X[n][0] for n in junction]))
    spec = {}
    for side in (-1, 1):
        J = [n for n in junction if (X[n][0] - x0) * side >= 0]
        Jx = np.array([X[n] for n in J])
        P = [p for p in rep if (X[p][0] - x0) * side > 0]
        L = min(np.linalg.norm(Jx - X[p], axis=1).min() for p in P)
        spec[side] = (J, float(L))
    return spec


if __name__ == '__main__':
    os.makedirs(SCRATCH, exist_ok=True)
    for nm in (NAME, PRE, MID):
        assert not os.path.exists(os.path.join(RUNS, nm)), f'{nm} exists; not overwriting'
    orig = os.path.join(RUNS, BASE, BASE + '.feb')
    spec = flare_spec(orig)
    print('flare: junction nodes / fade length per side', {k: (len(v[0]), round(v[1], 1)) for k, v in spec.items()})
    m0 = Model7(orig)
    rep = wall_reshape.reshape(m0, flare=spec)
    mid = os.path.join(SCRATCH, NAME + '_reshaped.feb')
    m0.write(mid)
    m = Model113(mid)
    m.log.append(f'{b135.LABEL}: wall_reshape (PVW flared toward the PeB over {[round(v[1], 1) for v in spec.values()]} mm, '
                 f'the PeB not narrowed): {rep}')
    b135.NAME = NAME
    b135.wrap2(m, orig)
    print(m.log[-1][:300])
    seg_up2(m)
    step2(INNER + OUTER_BOTTOM)(m)
    b135.walls_la(m)
    print(m.log[-1][:300])
    emit(PRE, m)
    run_main(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_batch129.py'), {'BASE': PRE, 'NAME': MID})
    run_main(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_batch129b.py'),
             {'SRC': MID, 'OLD': PRE, 'NAME': NAME})
    out = os.path.join(RUNS, NAME)
    notes = [open(os.path.join(RUNS, nm, nm + '.feb.changes.txt'), encoding='utf-8').read() for nm in (PRE, MID)]
    notes.append(open(os.path.join(out, NAME + '.feb.changes.txt'), encoding='utf-8').read())
    with open(os.path.join(out, NAME + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
        fo.write(f'base: {orig} (tools/build_batch136.py: reshaped and wrapped as {PRE}, then PM_conn on the inner arc as '
                 f'{MID}, then its law scaled)\n' + '\n'.join(notes))
    for nm in (PRE, MID):
        shutil.move(os.path.join(RUNS, nm), os.path.join(SCRATCH, nm))
    print('wrote', os.path.join(out, NAME + '.feb'))
