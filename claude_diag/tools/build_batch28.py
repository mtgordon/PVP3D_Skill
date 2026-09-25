"""Batch 28 (2026-09-24, user-approved to start at the beginning of the next session): the equilibrium check.
  L14_p3_hold  L12_stabdamp_p3 (L11_stab + mass damping C = 20/s from t = 0 + every surface pressure at 1/3; t = 1 with 0
               failed attempts in 11.5 min, LA 6.3 / 10.2 / 13.5 mm) run on to t = 2.0 with the load held (Amp-1 extends
               constant past t = 1; the damping stays on). Run length only. Why: the slow-ramp run L13_ramp3 showed the
               1 s dynamic runs lag the load (at 65 % load its LA was 2.5x the 1 s run's), so is the t = 1 state an
               equilibrium or still growing?
usage: py -3.10 build_batch28.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L12_stabdamp_p3', 'L12_stabdamp_p3.feb')
WANT = set(sys.argv[1:])

if not WANT or 'L14_p3_hold' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L14_p3_hold')), 'L14_p3_hold exists; not overwriting'
    m = Model6(BASE)
    amp = next(l for l in m.root.find('LoadData') if l.get('name') == 'Amp-1')
    assert amp.find('extend').text.upper() == 'CONSTANT', 'the load would not be held past t = 1'
    m.set_end_time(2.0)
    m.log.append('the load is held from t = 1 to 2 (Amp-1 extends CONSTANT); the damping (on from t = 0) stays on')
    emit('L14_p3_hold', m)
