"""Batch 39 (2026-09-25, user): batch 37-38 again with the LA pressure at 0.01 MPa "to get it more out of the way".

L27 (Load-LA 0.005, the other four 0.03 MPa) crawled from t ~ 0.32: the doubled wall pressure pushed the posterior wall
onto the LA early (PVW_LA touching at 8-13 facets by t 0.35). The user cancelled it and its four variants and asked for
Load-LA 0.01 MPa (10 kPa, 0.71x the source 0.014), which pushes the LA away from the PVW. NOT IN SOURCE (all five).
  L29_lofts_LA10kPa_rest30kPa  L26_lofts_la3_pm (the lofts line with PM_Plane) + Load-LA 0.01, the other four 0.03 MPa
and one change each on it, as batch 38:
  L29_lofts1030_refs10  max_refs 25 -> 10 (solver only)
  L29_lofts1030_pen01   PVW_LA penalty 0.5 -> 0.1 (auto_penalty kept; the contact pair stays)
  L29_lofts1030_opt25   opt_iter 10 -> 25 (solver only)
  L29_lofts1030_ramp3   the load ramped over 3 s (NOT IN SOURCE: the source ramps over exactly 1 s)
usage: py -3.10 build_batch39.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L26_lofts_la3_pm', 'L26_lofts_la3_pm.feb')
CONTROL = 'L29_lofts_LA10kPa_rest30kPa'
NEW = {'Load-LA': 0.01, 'Load-AVW': 0.03, 'Load-PVW': 0.03, 'Load-PeB-top': 0.03, 'Load-top': 0.03}
WANT = set(sys.argv[1:])


def pressures(m):
    done = set()
    for sl in m.root.find('Loads'):
        if sl.get('name') in NEW:
            p = sl.find('pressure')
            old = p.text
            p.text = '%g' % NEW[sl.get('name')]
            done.add(sl.get('name'))
            m.log.append(f'NOT IN SOURCE (user): {sl.get("name")} pressure {old} -> {p.text} MPa '
                         f'({NEW[sl.get("name")] / 0.014:.3g}x the source 0.014)')
    assert done == set(NEW), done


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
    m.log.append(f'NOT IN SOURCE: the load ramps over {s:g} s (Abaqus: 1 s); ends at t = '
                 f'{int(ctrl.find("time_steps").text) * float(ctrl.find("step_size").text):g}')


def pen01(m):
    assert m.set_contact('PVW_LA', penalty=0.1) == 1
    m.log.append('PVW_LA: penalty 0.5 -> 0.1 (auto_penalty 1 kept); the source contact pair is kept')


def build(name, base, fn):
    if WANT and name not in WANT:
        return
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(base)
    fn(m)
    emit(name, m)


build(CONTROL, BASE, pressures)
CBASE = os.path.join(RUNS, CONTROL, CONTROL + '.feb')
build('L29_lofts1030_refs10', CBASE, lambda m: m.set_solver(max_refs=10))
build('L29_lofts1030_pen01', CBASE, pen01)
build('L29_lofts1030_opt25', CBASE, lambda m: m.set_solver(opt_iter=25))
build('L29_lofts1030_ramp3', CBASE, lambda m: stretch_time(m, 3.0))
