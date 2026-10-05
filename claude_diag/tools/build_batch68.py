"""Batch 68 (2026-09-26 ~12:22): where does the best set settle? Hold the load past t = 1.

L56_springs_pc10_cu30_antside0 reached t = 1 (normal, 11 failed): +66.9 deg at t = 1, but +72.6 at t 0.98 and still
turning back at ~300 deg/s when the ramp ended (a dynamic run: the t = 1 state is not an equilibrium, batch 28). The
models already carry the settle phase (mass damping C = 20/s ramped in over t 1.0-1.05, loads held: curves extend
CONSTANT) but end at t = 1 (100 steps). One change (the run length; solver/run control, not physics):
  L61_springs_pc10_cu30_antside0_hold   L56_springs_pc10_cu30_antside0 run on to t = 1.5 (time_steps 100 -> 150)
Read with tools/hold_settle.py (--t0 1.0) and tools/peb_rotation.py --all.
usage: py -3.10 build_batch68.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def hold_to(m, t_end=1.5):
    ctrl = m.root.find('Control')
    dt = float(ctrl.find('step_size').text)
    el = ctrl.find('time_steps')
    old = el.text
    el.text = str(int(round(t_end / dt)))
    names = [b.get('name') for b in m.root.find('Loads').findall('body_load')]
    assert 'settle_damping' in names, names
    m.log.append(f'run control: time_steps {old} -> {el.text} (on to t = {t_end:g} with the loads held and the settle '
                 f'phase\'s mass damping, already in the model from t 1.0)')


BUILDS = (('L61_springs_pc10_cu30_antside0_hold', 'L56_springs_pc10_cu30_antside0'),)
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        hold_to(m)
        emit(name, m)
