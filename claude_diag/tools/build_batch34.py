"""Batch 34 (2026-09-25 morning, the user's choices after the overnight session): bigger time steps on the two lines.

The two lines the user wants tracked (both fast base, real cutback, PVW-LA penalty 0.5, no stab springs, no damping):
  springs  "the Abaqus replica" L19_pen5: every Abaqus CONN3D2 connector a FEBio spring, no lofts (closest to the source)
  lofts    the reference L21_lofts4_pen_r2: the user's four fitted lofts (AVW-Para, CL, USL, PM), the soft lofts as their
           connectors
Load-LA stays at 1/3 (0.014 -> 0.00467 MPa, NOT IN SOURCE: the user's choice until the LA's large displacement is
understood); the other four pressures stay at the source 0.014. Bigger time steps = dtmax 0.005 -> 0.01 or 0.02
(solver only; step_size 0.01 x 100 steps still ends at t = 1, the auto-stepper grows the step up to dtmax, so the floor
is 200 / 100 / 50 steps). One change per line from its _la3 baseline:
  L24_springs_la3        L19_pen5 + Load-LA x 1/3
  L24_springs_la3_dt01   + dtmax 0.01
  L24_springs_la3_dt02   + dtmax 0.02
  L24_lofts_la3          L21_lofts4_pen_r2 + Load-LA x 1/3
  L24_lofts_la3_dt01     + dtmax 0.01
  L24_lofts_la3_dt02     + dtmax 0.02
(The first build of this batch, with other names and one full-LA variant, is unrun in claude_diag/superseded/.)
usage: py -3.10 build_batch34.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASES = {'springs': os.path.join(RUNS, 'L19_pen5', 'L19_pen5.feb'),
         'lofts': os.path.join(RUNS, 'L21_lofts4_pen_r2', 'L21_lofts4_pen_r2.feb')}
WANT = set(sys.argv[1:])

for line, base in BASES.items():
    for dt in (None, 0.01, 0.02):
        name = f'L24_{line}_la3' + ('' if dt is None else f'_dt{int(round(dt * 1000)):02d}'.replace('_dt10', '_dt01')
                                     .replace('_dt20', '_dt02'))
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(base)
        m.log.append(f'the {line} line ({"the Abaqus replica" if line == "springs" else "the reference"})')
        m.scale_surface_load('Load-LA', 1 / 3)
        if dt is not None:
            m.set_solver(dtmax=dt)
            m.log.append(f'solver only: bigger time steps: the auto-stepper may take steps up to {dt:g} (was 0.005): at '
                         f'least {int(round(1 / dt))} steps to t = 1 instead of 200')
        emit(name, m)
