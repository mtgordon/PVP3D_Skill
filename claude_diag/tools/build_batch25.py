"""Batch 25 (2026-09-24, user-chosen), on the reference L9_fitall_edge (fast base: line search, every loft fitted, the
PM_PeB edge BC, arcus chain mass x10, Load-LA x 1/3), to t = 1.0, one change each:
  L11_rho   every loft at tissue density (1.06e-9): the lofts weigh 1355 g in L9_fitall_edge (AVW-Para, CL, USL and PM
            at the old pipe beam's 7.8e-07, P-arcus at 1.99e-07; the near-zero lofts already tissue), ~2.5 g at tissue
            density; the Abaqus connectors they replace are massless
  L11_stab  the 145 remaining zero-length stab_* ground springs removed (converted file, not in Abaqus; 0.02 N/mm each:
            56 on the perineal body (_PickedSet66), 42 on the apical vaginal wall (_PickedSet346), 47 on orphan nodes)
usage: py -3.10 build_batch25.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from variants2 import remove_stab_springs

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L9_fitall_edge', 'L9_fitall_edge.feb')
WANT = set(sys.argv[1:])


def new(name):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(BASE)


if not WANT or 'L11_rho' in WANT:
    m = new('L11_rho'); m.set_loft_density(1.06e-9); emit('L11_rho', m)
if not WANT or 'L11_stab' in WANT:
    m = new('L11_stab'); remove_stab_springs(m); emit('L11_stab', m)
