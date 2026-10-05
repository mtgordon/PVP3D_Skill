"""Batch 65 (2026-09-26 ~11:37): through the snap with a slower load ramp.

Past ~50 deg the runs snap and crawl: L56_springs_pc0_cu30 (P-arcus cut, CL/USL 30 %) jumps +39.9 -> +56.7 deg between
t 0.364 and 0.394 (peak +58.5), then the PVW runs at up to 460 mm/s and it crawls (32 failed at t 0.49);
L54_springs_pc10_cu30_LAfull peaks +53.1 at t 0.60, then the perineal body snaps at 600 mm/s. Abaqus Explicit rides
through a snap on inertia; an implicit run needs tiny steps. The 3 s ramp carried the soft-PVW runs further (batch 43).
One change each (the load curves and the step size x3, as build_batch43.ramp3; NOT IN SOURCE: the source ramps over 1 s):
  L58_springs_pc0_cu30_ramp3          L56_springs_pc0_cu30 + the 3 s ramp
  L58_springs_pc10_cu30_LAfull_ramp3  L54_springs_pc10_cu30_LAfull + the 3 s ramp
usage: py -3.10 build_batch65.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def ramp3(m, s=3.0):
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
        m.log.append(f'{path}: {old} -> {el.text}')
    m.log.append(f'NOT IN SOURCE: the load ramps over {s:g} s (Abaqus: 1 s); ends at t = '
                 f'{int(ctrl.find("time_steps").text) * float(ctrl.find("step_size").text):g}')


BUILDS = (('L58_springs_pc0_cu30_ramp3', 'L56_springs_pc0_cu30'),
          ('L58_springs_pc10_cu30_LAfull_ramp3', 'L54_springs_pc10_cu30_LAfull'))
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        ramp3(m)
        emit(name, m)
