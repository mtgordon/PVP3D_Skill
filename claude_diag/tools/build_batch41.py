"""Batch 41 (2026-09-25, user away ~1 h: "I want the PVW to go much further forward and the AVW to go further forward as
well. Work on figuring out how to make that happen by changing the material properties").

Forward = +y (anterior): Load-PVW pushes the posterior wall along (0, +0.83, +0.56), Load-AVW the anterior wall along
(0, -0.83, -0.56) (tools/load_direction.py). In the base L31_lofts_clusl30_sph50_LA10kPa (lofts line, CL/USL 30 %,
sphincter connectors 50 %, Load-LA 0.01, the other four 0.014) the mean displacement at t = 1 is PVW (0, -5.8, -20.4) and
AVW (0, +1.5, -22.8) mm: the PVW drifts backward as it descends. One change each, all NOT IN SOURCE:
  L32_pvw50     Vagina_PVW (2-term Ogden) c1, c2, k x 0.5: a softer posterior wall bulges forward under its pressure
  L32_pvw20     Vagina_PVW x 0.2
  L32_parcus50  the 26 P-arcus connectors (PVW to the posterior arcus chains) at 50 % force (spring scale 0.5)
  L32_peb20     the perineal body (PeB-Vagina500%stiffer, Ogden) c1, c2, k x 0.2 (to about the wall's stiffness)
usage: py -3.10 build_batch41.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L31_lofts_clusl30_sph50_LA10kPa', 'L31_lofts_clusl30_sph50_LA10kPa.feb')
WANT = set(sys.argv[1:])


def scale_ogden(m, name, f):
    mat = next(x for x in m.root.find('Material') if x.get('name') == name)
    assert mat.get('type') == 'Ogden', mat.get('type')
    changed = []
    for k in ('k', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6'):
        el = mat.find(k)
        if el is not None and float(el.text) != 0:
            old = el.text
            el.text = '%.6g' % (float(old) * f)
            changed.append(f'{k} {old} -> {el.text}')
    m.log.append(f'NOT IN SOURCE (material): {name} at {100 * f:g} % stiffness (Ogden c_i and k x {f:g}; m_i unchanged): '
                 + ', '.join(changed))


def scale_spring(m, name, f):
    dm = next(d for d in m.root.find('Discrete').findall('discrete_material') if d.get('name') == name)
    old = dm.find('scale').text
    dm.find('scale').text = '%g' % (float(old) * f)
    m.log.append(f'NOT IN SOURCE (connector law): {name} at {100 * f:g} % force: scale {old} -> {dm.find("scale").text}')


VARIANTS = {
    'L32_pvw50': lambda m: scale_ogden(m, 'Vagina_PVW', 0.5),
    'L32_pvw20': lambda m: scale_ogden(m, 'Vagina_PVW', 0.2),
    'L32_parcus50': lambda m: scale_spring(m, 'Parcus_conn_mat', 0.5),
    'L32_peb20': lambda m: scale_ogden(m, 'PeB-Vagina500%stiffer', 0.2),
}
for name, fn in VARIANTS.items():
    if WANT and name not in WANT:
        continue
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(BASE)
    fn(m)
    emit(name, m)
