"""Batch 21 (2026-09-24, task 3: the in-situ check of the fitted lofts; the near-zero lofts), on the new base L6_ls
(user, after batch 20: L5X_all + line search on, lstol 0.9; it reached t = 1 with 6 failed attempts vs 99 for the
control, same solution), ending at t = 1.0:
  L7_conn   every non-P-arcus loft replaced by the Abaqus connectors it stands in for, as FEBio nonlinear springs
            (variants6.lofts_to_connectors; the P-arcus lofts stay: fitted, and checked in situ in batch 17). The
            reference for L6_fitall: each family's connector-end elongations there should match these.
  L7_nzconn only the five near-zero lofts (PM_PeB x2, PM_avw_bottom x2, PeB-constrin) -> their connectors as springs:
            the faithful form of L6_pmpeb + L6_pebc (a near-zero loft is a soft membrane with hundreds of interior
            nodes; the connectors it stands in for are 1D and massless)
  L7_pebc   L6_pebc (PeB-constrin its own near-zero fit, tissue density; failed at t = 0.270 with lstol 0) with the
            line search on: does it get past that wall?
(First built on L5X_all with lstol 0, not run: claude_diag/superseded/L7_*_lstol0_unrun/.)
usage: py -3.10 build_batch21.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit, LOFT_FAMILIES
from loft_fit_all import fit_all

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L6_ls', 'L6_ls.feb')
TISSUE = 1.06e-9
WANT = set(sys.argv[1:])


def want(name):
    return not WANT or name in WANT


def new(name, base=BASE, t_end=None):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(base)
    if t_end:
        m.set_end_time(t_end)
    return m


if want('L7_conn'):
    m = new('L7_conn'); m.lofts_to_connectors(LOFT_FAMILIES); emit('L7_conn', m)
if want('L7_nzconn'):
    m = new('L7_nzconn')
    m.lofts_to_connectors(('PM_PeB_Left_', 'PM_PeB_Right_', 'PM_avw_bottom_left_', 'PM_avw_bottom_right_', 'PeB-constrin'))
    emit('L7_nzconn', m)
if want('L7_pebc'):
    r = fit_all(BASE, ['PeB-constrin_fan'])['PeB-constrin_fan']
    m = new('L7_pebc')
    m.set_loft_material('PeB-constrin_fan', r['c'], r['m'], 250.0, density=TISSUE,
                        note=f"loft_fit_all.py strip fit to its {r['n']} Abaqus connectors, rms {100 * r['err']:.1f} % "
                             f"over u = 2-47 mm")
    emit('L7_pebc', m)
