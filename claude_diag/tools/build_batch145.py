"""Batch 145+ (2026-10-02, the user: reproduce the paper's P1 / P2 (Luo et al. 2015) by changing the "healthy" connective
tissue and vaginal walls of the new base L144_seamspr_pm_allcontact). Builds, for one set of healthy changes, the healthy
model and the paper's cases P1 / P2 on it (build_batch107 P1 / P2: they MULTIPLY the healthy <scale>s, so the impairments
scale from the new healthy values; the walls get the same scale in all three).
Healthy changes (each NOT FAITHFUL; bounds x0.2 to x5):
  apical=S   every CL / USL spring set's <scale> x S (CL-L/R_conn_1-3, USL-L/R_conn_1-2)
  para=S     AVW-Para-L / R_conn x S                 parcus=S   Parcus_conn x S
  pmconn=S   PM_conn + PM-LA-x%stiff_mat x S          sph=S      LA_sphincter_side_conn + LA_sphincter_post_conn x S
  avw=S      Vagina_AVW + Vagina_Cervix Yeoh c1, c2, k x S
  pvw=S      Vagina_PVW Yeoh c1, c2, k x S
  walls=S    all three walls (AVW, cervix, PVW) x S: the user's choice for the fit (2026-10-02: "keep the walls the same,
             one shared scale"); not combined with avw / pvw
usage: py -3.10 build_batch145.py BATCH TAG [key=S ...] [--cases healthy,P1,P2] [--base RUN]
  -> runs/L<BATCH>_<TAG>_<case>/   (an existing run folder is never overwritten)"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch44 import spring_scale  # noqa: E402
from build_batch50 import CLUSL  # noqa: E402
from build_batch107 import P1, P2  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L144_seamspr_pm_allcontact'
LABEL = 'NOT FAITHFUL (the paper fit: healthy connective tissue / walls)'
FAMILIES = {'apical': CLUSL, 'para': ('AVW-Para-L_conn', 'AVW-Para-R_conn'), 'parcus': ('Parcus_conn',),
            'pmconn': ('PM_conn', 'PM-LA-x%stiff_mat'), 'sph': ('LA_sphincter_side_conn', 'LA_sphincter_post_conn')}
WALLS = {'avw': ('Vagina_AVW', 'Vagina_Cervix'), 'pvw': ('Vagina_PVW',),
         'walls': ('Vagina_AVW', 'Vagina_Cervix', 'Vagina_PVW')}   # the user, 2026-10-02: keep the walls the same (one scale)
CASES = {'healthy': None, 'P1': P1, 'P2': P2}


def wall_scale(m, mats, s):
    for nm in mats:
        mat = m._material(nm)
        old = {k: mat.find(k).text for k in ('c1', 'c2', 'k')}
        for k in ('c1', 'c2', 'k'):
            mat.find(k).text = '%.7g' % (float(old[k]) * s)
        m.log.append(f'{LABEL}: wall {nm} Yeoh c1, c2, k x {s:g}: {old} -> {({k: mat.find(k).text for k in old})}')


def healthy(m, scales):
    for key, s in scales.items():
        assert 0.2 - 1e-9 <= s <= 5 + 1e-9, f'{key}={s} outside x0.2-x5'
        if key in FAMILIES:
            for st in FAMILIES[key]:
                spring_scale(m, st, s, LABEL + f', {key}')
        else:
            wall_scale(m, WALLS[key], s)


def build(batch, tag, scales, cases=('healthy', 'P1', 'P2'), base=BASE):
    names = []
    for c in cases:
        name = f'L{batch}_{tag}_{c}'
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        m.log.append(f'{name}: {base} + healthy changes {scales or "none"} + case {c}')
        healthy(m, scales)
        if CASES[c]:
            CASES[c](m)
        emit(name, m)
        names.append(name)
    return names


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('batch')
    ap.add_argument('tag')
    ap.add_argument('kv', nargs='*')
    ap.add_argument('--cases', default='healthy,P1,P2')
    ap.add_argument('--base', default=BASE)
    a = ap.parse_args()
    sc = {k: float(v) for k, v in (x.split('=') for x in a.kv)}
    for k in sc:
        assert k in FAMILIES or k in WALLS, k
    assert not ('walls' in sc and ('avw' in sc or 'pvw' in sc)), 'walls is not combined with avw / pvw'
    print(' '.join(build(a.batch, a.tag, sc, a.cases.split(','), a.base)))
