"""Batch 29 (2026-09-25, overnight; the user: "keep working on getting the model to run well (stable, modeling it well,
and then fast)"): the stab-spring-free model at the SOURCE pressures, fast base, to t = 1.0.

Base "src" = L11_stab (the reference L9_fitall_edge without the 145 stab_* ground springs, which Abaqus does not have)
with Load-LA back at the source 0.014 MPa, so all five pressures are the Abaqus values again. Without the springs the
apex (USL lofts, _PickedSet64/346/347) moves at 100-200 mm/s and the 1/3-LA-load runs crawled near t = 0.5-0.6; with
every pressure at 1/3 plus mass damping the run was clean (L12_stabdamp_p3). One change each on src:
  L15_src        control
  L15_psym       the five pressure loads with their exact tangent: surface_load symmetric_stiffness 1 -> 0. The global
                 matrix is already unsymmetric (the sliding-elastic contacts ask for it: log "matrix format:
                 unsymmetric"), so this costs nothing per solve; only the pressures' follower stiffness was symmetrised.
                 Numerical only, the same physics
  L15_avwconn    the AVW-Para L/R lofts -> their 70 Abaqus connectors as springs (tension only, as in the source). In
                 situ the lofts are compressed (connector ends -1.2 to -1.7 mm at t 0.55-0.63), which the source's
                 connectors cannot carry; a thin shell in compression can wrinkle
  L15_cluslconn  the CL/USL L/R lofts -> their 82 connectors: the apex supports, where the fastest nodes are
  L15_conn       all 12 non-P-arcus lofts -> their 175 connectors (L7_conn's change; the most faithful supports)
  L15_rho        every loft at tissue density (1.06e-9; the lofts now weigh 1.36 kg against ~70 g of tissue, the
                 connectors they replace are massless). Faithful mass for the like-for-like dynamics
  L15_damp       mass damping (C = 20/s) on from t = 0. NOT IN SOURCE (numerical)
  L15_damp_r3    L15_damp + the load ramped over 3 s (step_size, dtmax x3). NOT IN SOURCE: the user's "slow ramp,
                 damping, or both"; the undamped slow ramp L13_ramp3 snapped through at 81 % of the 1/3-LA load
Every run also plots velocity (output only).
usage: py -3.10 build_batch29.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L11_stab', 'L11_stab.feb')
PRESSURES = ('Load-AVW', 'Load-PVW', 'Load-PeB-top', 'Load-top', 'Load-LA')
AVW = ('AVW-Para-L', 'AVW-Para-R')
CLUSL = ('CL-L', 'CL-R', 'USL-L', 'USL-R')
WANT = set(sys.argv[1:])


def want(name):
    return not WANT or name in WANT


def src(name):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(BASE)
    for sl in m.root.find('Loads'):
        if sl.get('name') in PRESSURES:
            p = sl.find('pressure')
            if abs(float(p.text) - 0.014) > 1e-9:
                old = p.text
                p.text = '0.014'
                m.log.append(f'{sl.get("name")} pressure {old} -> 0.014 MPa: back at the Abaqus value (all five '
                             f'pressures are now the source loads)')
    assert all(abs(float(s.find('pressure').text) - 0.014) < 1e-9 for s in m.root.find('Loads')
               if s.get('name') in PRESSURES)
    plot = m.root.find('Output').find('plotfile')
    if not any(v.get('type') == 'velocity' for v in plot):
        ET.SubElement(plot, 'var', type='velocity')
        m.log.append('output only: velocity added to the plot file (no model change)')
    return m


def exact_pressure_tangent(m):
    for sl in m.root.find('Loads'):
        if sl.get('name') in PRESSURES:
            s = sl.find('symmetric_stiffness')
            old = s.text
            s.text = '0'
            m.log.append(f'{sl.get("name")}: symmetric_stiffness {old} -> 0 (the exact, unsymmetric follower-pressure '
                         f'tangent; the global matrix is already unsymmetric for the contacts). Numerical only')


def damping_from_start(m):
    lc = next(l for l in m.root.find('LoadData') if l.get('name') == 'settle_on')
    users = [e for e in m.root.iter() if e.get('lc') == lc.get('id')]
    assert [u.tag for u in users] == ['C'], [u.tag for u in users]
    pts = lc.find('points')
    old = [p.text for p in pts]
    for p in list(pts):
        pts.remove(p)
    for txt in ('0,1', '1.05,1'):
        ET.SubElement(pts, 'pt').text = txt
    m.log.append(f'NOT IN SOURCE (numerical): mass damping settle_damping (C = 20/s) on from t = 0: load curve settle_on '
                 f'{old} -> [0,1; 1.05,1] (it drives only the damping); Abaqus explicit has no such damping')


def stretch_time(m, s=3.0):
    ld = m.root.find('LoadData')
    for name in ('Amp-1', 'settle_on'):
        lc = next(l for l in ld if l.get('name') == name)
        pts = lc.find('points')
        old = [p.text for p in pts]
        for p in pts:
            t, v = p.text.split(',')
            p.text = f'{float(t) * s:g},{v}'
        m.log.append(f'load curve {name}: time x{s:g}: {old} -> {[p.text for p in pts]}')
    ctrl = m.root.find('Control')
    for path in ('step_size', 'time_stepper/dtmax'):
        el = ctrl.find(path)
        old = el.text
        el.text = '%g' % (float(old) * s)
        m.log.append(f'{path}: {old} -> {el.text} (the same number of steps over the {s:g}x longer ramp)')
    m.log.append(f'NOT IN SOURCE: the load ramps over {s:g} s (Abaqus: 1 s); time_steps {ctrl.find("time_steps").text} '
                 f'x step {ctrl.find("step_size").text} ends at t = '
                 f'{int(ctrl.find("time_steps").text) * float(ctrl.find("step_size").text):g}')


if want('L15_src'):
    m = src('L15_src'); m.log.append('control: the source pressures, no damping, 1 s ramp'); emit('L15_src', m)
if want('L15_psym'):
    m = src('L15_psym'); exact_pressure_tangent(m); emit('L15_psym', m)
if want('L15_avwconn'):
    m = src('L15_avwconn'); m.lofts_to_connectors(AVW); emit('L15_avwconn', m)
if want('L15_cluslconn'):
    m = src('L15_cluslconn'); m.lofts_to_connectors(CLUSL); emit('L15_cluslconn', m)
if want('L15_conn'):
    m = src('L15_conn'); m.lofts_to_connectors(); emit('L15_conn', m)
if want('L15_rho'):
    m = src('L15_rho'); m.set_loft_density(1.06e-9); emit('L15_rho', m)
if want('L15_damp'):
    m = src('L15_damp'); damping_from_start(m); emit('L15_damp', m)
if want('L15_damp_r3'):
    m = src('L15_damp_r3'); damping_from_start(m); stretch_time(m, 3.0); emit('L15_damp_r3', m)
