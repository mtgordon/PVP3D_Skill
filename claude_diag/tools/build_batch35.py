"""Batch 35 (2026-09-25 morning, user-approved): smaller steps and a stricter stepper on both lines, Load-LA x 1/3.

Batch 34 showed bigger time steps give the same solution, only slower (dtmax 0.01: 41 min vs 28 min on the lofts line;
0.02 behind everywhere): each big step needs more Newton iterations and fails more often. So the other direction, one
change each from the batch-34 baselines (L24_springs_la3: 40 min, 278 steps, 6868 iterations; L24_lofts_la3: 28 min,
215 steps, 2788 iterations; 5 threads each):
  L25_<line>_la3_dt0025  dtmax 0.005 -> 0.0025 (step_size 0.01 -> 0.0025 with time_steps 100 -> 400, still t = 1)
  L25_<line>_la3_opt10   opt_iter 25 -> 10 (the auto-stepper grows the step only while a step needs fewer than 10
                         iterations, and shrinks it above)
Solver only. <line> = springs (the Abaqus replica line) or lofts (the reference line).
usage: py -3.10 build_batch35.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

for line in ('springs', 'lofts'):
    base = os.path.join(RUNS, f'L24_{line}_la3', f'L24_{line}_la3.feb')
    for tag in ('dt0025', 'opt10'):
        name = f'L25_{line}_la3_{tag}'
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(base)
        if tag == 'dt0025':
            ctrl = m.root.find('Control')
            for k, v in (('step_size', '0.0025'), ('time_steps', '400')):
                old = ctrl.find(k).text
                ctrl.find(k).text = v
                m.log.append(f'{k} {old} -> {v} (still ends at t = 1)')
            m.set_solver(dtmax=0.0025)
            m.log.append('solver only: smaller time steps (dtmax 0.005 -> 0.0025)')
        else:
            m.set_solver(opt_iter=10)
            m.log.append('solver only: the auto-stepper targets 10 iterations per step instead of 25')
        emit(name, m)
