"""Batch 23 (2026-09-24, user: "for PM-PeB, make the boundary condition on it be all of the nodes along the edge"):
  L9_nzw_edge  L8_nzweak (the five near-zero lofts fitted, tissue density, line search on) + the PM_PeB lofts' whole
               anchor edge under BC-RP-PM-PeB-origins (variants6.fix_anchor_edge: + 7 left / 9 right free midpoint
               nodes between the 8 fixed origins; with the weak loft they swung up to 3.1-3.6 mm in L6_pmpeb)
usage: py -3.10 build_batch23.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def new(name, base):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(base)


if not WANT or 'L9_nzw_edge' in WANT:
    m = new('L9_nzw_edge', os.path.join(RUNS, 'L8_nzweak', 'L8_nzweak.feb'))
    for loft in ('PM_PeB_Left_fan', 'PM_PeB_Right_fan'):
        m.fix_anchor_edge(loft, 'BC-RP-PM-PeB-origins')
    emit('L9_nzw_edge', m)
