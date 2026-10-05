"""Batch 89 (2026-09-27 ~21:40, the user away ~8 h): push the faithful walls on from the batch-88 leader.

L85_newline_vwyeoh_rhoi0 (the new springs line + the Yeoh walls + solver rhoi 0) led batch 88 (t 0.726 at 21:17 vs the
base 0.647). Where it crawls (tools/contact_onset.py, t 0.70-0.73): steps of dt 1e-4..1e-3 with 10-50 iterations; the
failed and slow iterations meet the energy norm by ~5 orders and never the displacement norm (current 2e-4..5e-3 against
2.7e-6 required); the PVW's lateral edge row on the LA (1894-1898, 1984-1987, 5400/5401, 7331; 8-14 PVW_LA facets
touching at 0.005-0.015 MPa) reverses direction step to step at 0.2-0.8 m/s, and the canal seam (AVW 2076-2078, 2155-2158;
PVW 2571, 2588/2589, 7033) at 0.2-1 m/s. One change each from L85_newline_vwyeoh_rhoi0:
  L86_newline_vwyeoh_rhoi0_dtol01     solver dtol 0.001 -> 0.01 (the displacement norm is the only one not met; a
                                      tolerance, not physics: check the solution against the parent at matching t)
  L86_newline_vwyeoh_rhoi0_newton     solver max_ups 10 -> 0 (full Newton: the stiffness reformed every iteration, so it
                                      follows the contact state as it switches)
  L86_newline_vwyeoh_rhoi0_pvwpen025  PVW_LA penalty 0.5 -> 0.25 (5 -> 0.5 was the best fix so far for the PVW edge)
  L86_newline_vwyeoh_rhoi0_se1segup5  SlidingElastic1 seg_up 0 -> 5 (facets re-chosen only in each step's first 5
                                      iterations; seg_up 2 froze them too early on the canal, L75_nopinch_vwyeoh_se1segup2)
usage: py -3.10 build_batch89.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch79 import se1

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def ctrl(note, **kw):
    def fn(m):
        m.set_control(**kw)
        m.log.append(note)
    return fn


def pvw(**kw):
    def fn(m):
        assert m.set_contact('PVW_LA', **kw) == 1
    return fn


BUILDS = (('L86_newline_vwyeoh_rhoi0_dtol01',
           ctrl('solver only: dtol 0.001 -> 0.01 (convergence tolerance); the model is unchanged', dtol=0.01)),
          ('L86_newline_vwyeoh_rhoi0_newton',
           ctrl('solver only: max_ups 10 -> 0 (full Newton); the model is unchanged', max_ups=0)),
          ('L86_newline_vwyeoh_rhoi0_pvwpen025', pvw(penalty=0.25)),
          ('L86_newline_vwyeoh_rhoi0_se1segup5', se1(seg_up=5)))
# ~22:35: dtol 0.01 got furthest (t 0.808 at 22:30, the solution within 0.2 mm of the parent's at t 0.76); canal seg_up 5
# the next best setting (steady, 2 failed). The two together, one change from L86_newline_vwyeoh_rhoi0_dtol01:
#   L88_newline_vwyeoh_rhoi0_dtol01_se1segup5
LATER = (('L88_newline_vwyeoh_rhoi0_dtol01_se1segup5', 'L86_newline_vwyeoh_rhoi0_dtol01', se1(seg_up=5)),)
if __name__ == '__main__':
    base = 'L85_newline_vwyeoh_rhoi0'
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
    for name, b, fn in LATER:
        if name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, b, b + '.feb'))
        fn(m)
        emit(name, m)
