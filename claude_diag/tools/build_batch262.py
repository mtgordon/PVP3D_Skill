"""Batch 262 (2026-10-10, the user: "go ahead and do a loft + springs for the tube as well to see how it does"). The springs +
lofts version of the tube main (L258_tubecand_*; memory feedback-three-base-versions, version 2): every connector family
as a chain of springs inside a weak sheet (build_batch111 chains_all, sheet 5 % of the family's fibre E1, as
L111_chains_all_m05 / L132_chains_tube_pm), built from the tube main's HEALTHY file; a case's impairment scales both the
chain springs and the family's sheet (as L111_chains_all_m05_paperP1), the levator by la_healthy_pct.
The tube's curve_end_PM_conn springs share PM_conn_mat; they get their own copy first (curve_end_PM_conn_mat, same table),
so PM_conn can become chains (the curve-end springs stay single springs, as the other PM-arch ones are kept by their rule).
Cases (build_batch244): healthy, MAX (levator 40 %, CL / USL x0.2, AVW-Para x0.5, Parcus x0.05).
NOT IN SOURCE (the chains' interior nodes and the sheets).
2026-10-10 (the user: "for ledge, don't add compression on the springs yet ... You can try it with the loft + springs"):
--src / --prefix build the same from another healthy file (e.g. the ledge's L243_ledge_cx3_src_healthy).
--skip FAM[,FAM]: families kept as single springs (e.g. AVW-Para).
--sheet M: the sheets' modulus as a fraction of each family's fibre E1 (default 0.05).
usage: py -3.10 build_batch262.py CASE [CASE ...] [--src HEALTHY_RUN --prefix NAME_PREFIX --sheet M]"""
import copy
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch111 import Model111, SHEETS  # noqa: E402
from build_batch107 import la_healthy_pct  # noqa: E402
from build_batch244 import CASES  # noqa: E402
from build_batch44 import spring_scale  # noqa: E402
from variants2 import renumber_discrete_materials  # noqa: E402
from variants7 import JOBS, emit  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L258_tubecand_healthy'


def own_curve_end_mat(m):
    disc = m.root.find('Discrete')
    mats = disc.findall('discrete_material')
    d = next(e for e in disc.findall('discrete') if e.get('discrete_set') == 'curve_end_PM_conn')
    src = mats[int(d.get('dmat')) - 1]
    assert src.get('name') == 'PM_conn_mat', src.get('name')
    new = copy.deepcopy(src)
    new.set('name', 'curve_end_PM_conn_mat')
    new.set('id', str(len(mats) + 1))
    disc.insert(list(disc).index(mats[-1]) + 1, new)
    d.set('dmat', str(len(mats) + 1))
    m.log.append('curve_end_PM_conn: its own copy of PM_conn_mat (curve_end_PM_conn_mat, same table) so PM_conn can be chained')


def build(case, src=SRC, prefix='L262_tubemain_chains', sheet_m=0.05, skip=()):
    la, ap, ant, post = CASES[case]
    m = Model111(os.path.join(RUNS, src, src + '.feb'))
    m.log.append(f'{src} -> the springs + lofts (chains) version, case {case}')
    if any(e.get('discrete_set') == 'curve_end_PM_conn' for e in m.root.find('Discrete').findall('discrete')):
        own_curve_end_mat(m)
    if la != 100:
        la_healthy_pct(m, la)
    key = {'AVW-Para': ant, 'CL': ap, 'USL': ap, 'P-arcus': post}
    gone = set()
    for sheet in SHEETS:
        s = next((v for k, v in key.items() if sheet.startswith(k)), 1.0)
        if any(sheet.startswith(k) for k in skip):
            for st in SHEETS[sheet][1]:
                if s != 1:
                    spring_scale(m, st, s, f'the case {case}')
            m.log.append(f'{sheet}: kept as single springs (--skip), case scale x{s:g}')
            continue
        gone.update(m.chains(sheet, sheet_m, s))
    disc = m.root.find('Discrete')
    for s in sorted(gone):
        pairs, mat, b, ds = m._spring_set(s)
        disc.remove(b)
        disc.remove(mat)
        m.mesh.remove(ds)
        m.log.append(f'{s}: its {len(pairs)} single springs removed (now chains)')
    renumber_discrete_materials(m)
    name = f'{prefix}_{case}'
    assert not os.path.exists(os.path.join(RUNS, name)), name
    emit(name, m)
    return name


if __name__ == '__main__':
    import argparse
    ap_ = argparse.ArgumentParser()
    ap_.add_argument('cases', nargs='+'); ap_.add_argument('--src', default=SRC); ap_.add_argument('--prefix', default='L262_tubemain_chains'); ap_.add_argument('--sheet', type=float, default=0.05); ap_.add_argument('--skip', default='')
    a = ap_.parse_args()
    for c in a.cases:
        print(build(c, a.src, a.prefix, a.sheet, tuple(x for x in a.skip.split(',') if x)))
