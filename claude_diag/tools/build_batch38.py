"""Batch 38 (2026-09-25, user-chosen): getting L27_lofts_LA5kPa_rest30kPa (the lofts line with Load-LA 0.005 MPa and
the other four pressures 0.03 MPa, NOT IN SOURCE) through its early PVW-LA contact.

L27 crawled from t ~ 0.32 (steps ~5e-5): the doubled wall pressure pushes the posterior wall onto the LA early (PVW_LA
touching at 8-13 facets, p up to 0.013 MPa by t 0.35; at the 1/3 loads it first touched after t ~ 0.95); converged steps
need 7-47 Newton iterations, and every failed attempt spends the whole 250 (max_refs 25 x max_ups 10). One change each:
  L28_lofts530_refs10  max_refs 25 -> 10: a failed attempt costs ~100 iterations before the step is halved. Solver only
  L28_lofts530_pen01   PVW_LA penalty 0.5 -> 0.1 (auto_penalty kept; the contact pair stays)
  L28_lofts530_opt25   opt_iter 10 -> 25 (back to the default): bigger steps allowed when each needs > 10 iterations
  L28_lofts530_ramp3   the load ramped over 3 s instead of 1 s (Amp-1 x3; step_size, dtmax x3; the settle damping,
                       off during the ramp, switches on at 3.0-3.15). NOT IN SOURCE: the source ramps over exactly 1 s
usage: py -3.10 build_batch38.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L27_lofts_LA5kPa_rest30kPa', 'L27_lofts_LA5kPa_rest30kPa.feb')
WANT = set(sys.argv[1:])


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


def variant(name, fn):
    if WANT and name not in WANT:
        return
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(BASE)
    fn(m)
    emit(name, m)


def pen01(m):
    assert m.set_contact('PVW_LA', penalty=0.1) == 1
    m.log.append('PVW_LA: penalty 0.5 -> 0.1 (auto_penalty 1 kept); the source contact pair is kept')


variant('L28_lofts530_refs10', lambda m: m.set_solver(max_refs=10))
variant('L28_lofts530_pen01', pen01)
variant('L28_lofts530_opt25', lambda m: m.set_solver(opt_iter=25))
variant('L28_lofts530_ramp3', lambda m: stretch_time(m, 3.0))
