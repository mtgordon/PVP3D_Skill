"""Batch 40 (2026-09-25, user): "make the CL and USL 50 % of their previous strength ... to get more movement in the
vagina walls", on the lofts line, with the vaginal-wall / cervix / perineal-body pressures back at the source 0.014 MPa
and Load-LA 0.01 MPa (user). Two changes from L26_lofts_la3_pm (the lofts line with PM_Plane; Load-LA 0.00467), both
NOT IN SOURCE:
  - the CL-L, CL-R, USL-L, USL-R lofts at 50 % strength: their fitted 1-term Ogden c1 x 0.5 and k x 0.5 (k stays 250 mu0),
    which halves the loft stress at every stretch (Ogden stress is linear in c1);
  - Load-LA 0.00467 -> 0.01 MPa (0.71x the source 0.014); the other four stay at 0.014.
  L30_lofts_clusl50_LA10kPa
usage: py -3.10 build_batch40.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L26_lofts_la3_pm', 'L26_lofts_la3_pm.feb')
LOFTS = ('CL-L_fan', 'CL-R_fan', 'USL-L_fan', 'USL-R_fan')
WANT = set(sys.argv[1:])

name = 'L30_lofts_clusl50_LA10kPa'
if not WANT or name in WANT:
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(BASE)
    doms = {d.get('name'): d.get('mat') for d in m.root.find('MeshDomains')}
    mats = {x.get('name'): x for x in m.root.find('Material')}
    for loft in LOFTS:
        mat = mats[doms[loft]]
        assert mat.get('type') == 'Ogden', mat.get('type')
        assert all(float(mat.find(f'c{i}').text) == 0 for i in range(2, 7) if mat.find(f'c{i}') is not None)
        old = {k: mat.find(k).text for k in ('c1', 'k')}
        for k in ('c1', 'k'):
            mat.find(k).text = '%.6g' % (float(old[k]) * 0.5)
        m.log.append(f'NOT IN SOURCE (user): {loft} loft ({mat.get("name")}) at 50 % strength: c1 {old["c1"]} -> '
                     f'{mat.find("c1").text}, k {old["k"]} -> {mat.find("k").text} (k stays 250 mu0; m1 unchanged)')
    sl = next(s for s in m.root.find('Loads') if s.get('name') == 'Load-LA')
    old = sl.find('pressure').text
    sl.find('pressure').text = '0.01'
    m.log.append(f'NOT IN SOURCE (user): Load-LA pressure {old} -> 0.01 MPa (0.71x the source 0.014); the other four '
                 f'pressures stay at the source 0.014')
    emit(name, m)

# The matching control (user-approved): the same pressures, CL/USL at full strength: one change from L26_lofts_la3_pm.
ctrl = 'L30_lofts_LA10kPa'
if not WANT or ctrl in WANT:
    assert not os.path.exists(os.path.join(RUNS, ctrl)), f'{ctrl} exists; not overwriting'
    m = Model6(BASE)
    sl = next(s for s in m.root.find('Loads') if s.get('name') == 'Load-LA')
    old = sl.find('pressure').text
    sl.find('pressure').text = '0.01'
    m.log.append(f'NOT IN SOURCE (user): Load-LA pressure {old} -> 0.01 MPa (0.71x the source 0.014); the other four '
                 f'pressures stay at the source 0.014. The control of L30_lofts_clusl50_LA10kPa (CL/USL at full strength)')
    emit(ctrl, m)

# L30_lofts_clusl30_LA10kPa (user, 2026-09-25): CL/USL at 30 % of the original fitted strength (c1 and k x 0.3), the
# same pressures (Load-LA 0.01, the rest 0.014). Two changes from L26_lofts_la3_pm, one from the control L30_lofts_LA10kPa.
n30 = 'L30_lofts_clusl30_LA10kPa'
if not WANT or n30 in WANT:
    assert not os.path.exists(os.path.join(RUNS, n30)), f'{n30} exists; not overwriting'
    m = Model6(BASE)
    doms = {d.get('name'): d.get('mat') for d in m.root.find('MeshDomains')}
    mats = {x.get('name'): x for x in m.root.find('Material')}
    for loft in LOFTS:
        mat = mats[doms[loft]]
        old = {k: mat.find(k).text for k in ('c1', 'k')}
        for k in ('c1', 'k'):
            mat.find(k).text = '%.6g' % (float(old[k]) * 0.3)
        m.log.append(f'NOT IN SOURCE (user): {loft} loft ({mat.get("name")}) at 30 % strength: c1 {old["c1"]} -> '
                     f'{mat.find("c1").text}, k {old["k"]} -> {mat.find("k").text} (k stays 250 mu0; m1 unchanged)')
    sl = next(s for s in m.root.find('Loads') if s.get('name') == 'Load-LA')
    old = sl.find('pressure').text
    sl.find('pressure').text = '0.01'
    m.log.append(f'NOT IN SOURCE (user): Load-LA pressure {old} -> 0.01 MPa (0.71x the source 0.014); the other four '
                 f'pressures stay at the source 0.014')
    emit(n30, m)

# L31_lofts_clusl30_sph50_LA10kPa (user, 2026-09-25): "do the 30 % of CL/USL and 50 % on all of the LA_PCMPRM" = the 40
# sphincter connectors from LA_PCMPRM to the perineal body (LA_sphincter_side_conn 16, LA_sphincter_post_conn 24; Abaqus
# CONN3D2) at 50 % force (nonlinear spring scale 1 -> 0.5), in one run with CL/USL at 30 % (user: one run with both).
# They carry the perineal body (4.2 N up, 2.7 N back at t = 1 in L30_lofts_clusl50_LA10kPa). NOT IN SOURCE.
n31 = 'L31_lofts_clusl30_sph50_LA10kPa'
if not WANT or n31 in WANT:
    assert not os.path.exists(os.path.join(RUNS, n31)), f'{n31} exists; not overwriting'
    m = Model6(os.path.join(RUNS, 'L30_lofts_clusl30_LA10kPa', 'L30_lofts_clusl30_LA10kPa.feb'))
    for dm in m.root.find('Discrete').findall('discrete_material'):
        if dm.get('name') in ('LA_sphincter_side_conn_mat', 'LA_sphincter_post_conn_mat'):
            old = dm.find('scale').text
            dm.find('scale').text = '0.5'
            m.log.append(f'NOT IN SOURCE (user): {dm.get("name")} (the sphincter connectors LA_PCMPRM - perineal body) '
                         f'at 50 % force: scale {old} -> 0.5')
    assert sum('at 50 % force' in l for l in m.log) == 2
    emit(n31, m)
