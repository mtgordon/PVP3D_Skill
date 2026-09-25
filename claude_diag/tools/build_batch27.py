"""Batch 27 (2026-09-24, the user's idea): a slower load ramp instead of damping, on L11_stab (= the reference
L9_fitall_edge without the stab_* ground springs), fast base:
  L13_ramp3  the load ramped over 3 s instead of 1 s: Amp-1 (all five pressures) smooth step 0 -> 1 over t = 0-3; the
             settle damping's switch-on moved from t = 1.0-1.05 to 3.0-3.15 (still off throughout the ramp, as in
             L11_stab); step_size 0.01 -> 0.03 and dtmax 0.005 -> 0.015 (the same number of steps), 100 steps to t = 3.
             Without the springs the posterior vaginal wall and the USL lofts move at 100-200 mm/s and the runs crawl
             near t = 0.5; a 3x slower ramp cuts the speeds ~3x and the inertial forces ~9x. NOT IN SOURCE: Abaqus ramps
             over exactly 1 s (*Dynamic, Explicit, time 1.0, Amp-1 SMOOTH STEP), so this is a route to the loaded state,
             not the like-for-like timing. Compare with the 1 s runs at the same load: t(ramp3) = 3 x t(1 s).
usage: py -3.10 build_batch27.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L11_stab', 'L11_stab.feb')
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
    m.log.append(f'NOT IN SOURCE: the load ramps over {s:g} s (Abaqus: 1 s); time_steps {ctrl.find("time_steps").text} '
                 f'x step {ctrl.find("step_size").text} ends at t = {int(ctrl.find("time_steps").text) * float(ctrl.find("step_size").text):g}')


if not WANT or 'L13_ramp3' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L13_ramp3')), 'L13_ramp3 exists; not overwriting'
    m = Model6(BASE); stretch_time(m, 3.0); emit('L13_ramp3', m)
