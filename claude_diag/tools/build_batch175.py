"""Batch 175 (2026-10-05, the user's 8-10 h unattended plan: reproduce the paper's P1 / P2 on the non-tube seam-band base by
tightening the healthy support tissues; no whole-wall stiffening, no new geometry).
The base L175_band3_pmd_sph025 = L149_seamspr_pm_outer with:
  1. the canal seam springs removed (constraint canal_seam_springs; the seam is contact-only, as in the source);
  2. the 3 mm seam band, Yeoh c1 x3 (build_batch100.band(3, 3); seam only, ~10 % of wall elements; NOT IN SOURCE);
  3. wall damping C 200/s on the AVW, cervix and PVW (build_batch165.wall_damping; NOT IN SOURCE);
  4. PM shell damping C 200/s (body_load "mass damping" elem_set PM_Plane, as L170; NOT IN SOURCE);
  5. LA_sphincter_side_conn_mat scale x0.25 (build_batch174.soften with FACTOR 0.25; the user's standing rule; NOT IN SOURCE).
Healthy base only (not run on its own until round 2); the fit cases come from
  py -3.10 build_batch145.py BATCH TAG key=S ... --cases P1 --base L175_band3_pmd_sph025
usage: py -3.10 build_batch175.py [BAND_MM [C1_FACTOR]]   (default 3 3; e.g. 6 -> L175_band6_pmd_sph025; 3 5 -> L175_band3x5_pmd_sph025)"""
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch100 import band  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch165 import wall_damping, C  # noqa: E402
import build_batch174  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L149_seamspr_pm_outer'
TMP = 'L175_band3_pmd_sph025_tmp'
NAME = 'L175_band3_pmd_sph025'
WALLS = ('_PickedSet347', '_PickedSet346', '_PickedSet64')


def pm_damping(t):
    assert t.count('<ShellDomain name="PM_Plane"') == 1
    add = f'\t\t<body_load name="pm_damping" type="mass damping" elem_set="PM_Plane">\n\t\t\t<C>{C}</C>\n\t\t</body_load>\n'
    i = t.index('\t</Loads>')
    return t[:i] + add + t[i:], (f'NOT IN SOURCE: mass damping on the PM shell (PM_Plane), C {C}/s from t = 0, as the walls\' '
                                 f'(as L170_avwcx5_pvwband3_pmdamp_P1)')


def make(name, bands, tag='tools/build_batch175.py'):
    """bands: [(width, factor, walls), ...] (build_batch100.band; each wall in one band only)."""
    tmp = name + '_tmp'
    for nm in (tmp, name):
        assert not os.path.exists(os.path.join(RUNS, nm)), nm
    m = Model7(os.path.join(RUNS, SRC, SRC + '.feb'))
    cons = m.root.find('Constraints')
    seam = [c for c in cons if c.get('name') == 'canal_seam_springs']
    assert len(seam) == 1
    n = len(seam[0].findall('linear_constraint'))
    cons.remove(seam[0])
    m.log.append(f'as the source: the canal seam springs removed (constraint canal_seam_springs, {n} linear constraints); the '
                 f'seam is contact-only (SlidingElastic1)')
    for w, f, walls in bands:
        band(w, f, walls)(m)
    emit(tmp, m)
    p = os.path.join(RUNS, tmp, tmp + '.feb')
    t = open(p, encoding='utf-8').read()
    msgs = []
    t, msg = wall_damping(t, WALLS)
    msgs.append(msg)
    t, msg = pm_damping(t)
    msgs.append(msg)
    build_batch174.FACTOR = 0.25
    t, msg = build_batch174.soften(t)
    msgs.append(msg)
    notes = open(p + '.changes.txt', encoding='utf-8').read()
    write(name, t, notes + '\n'.join(msgs) + f'\n({tag}; base {SRC})\n')
    shutil.rmtree(os.path.join(RUNS, tmp))


ALL = ('Vagina_AVW', 'Vagina_Cervix', 'Vagina_PVW')

if __name__ == '__main__':
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    F = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if (W, F) != (3, 3):
        NAME = f'L175_band{W}_pmd_sph025' if F == 3 else f'L175_band{W}x{F}_pmd_sph025'
    make(NAME, [(W, F, ALL)])
