"""Batch 22 (2026-09-24, user-approved after batches 20-21), all with the line search on, ending at t = 1.0.
Fast base L6_ls (L5X_all + lstol 0.9):
  L8_fitall   every non-P-arcus loft its own fitted material (as L6_fitall, which stalled at t = 0.908 without the
              line search): the loft counterpart of L7_conn (the connectors as springs) at t = 1
  L8_nzweak   the five near-zero lofts (PM_PeB x2, PM_avw_bottom x2, PeB-constrin) their own fits, tissue density:
              the loft counterpart of L7_nzconn
  L8_aggr     a real cutback: aggressiveness 0 -> 1 with cutback 0.5 (with aggressiveness 0, FEBio takes dt0/21 off the
              step at each retry and ignores cutback: FECore FETimeStepController::Retry)
Like-for-like base L6D_x (LPFmsTk_DM1 = Abaqus chain mass, + Load-LA x 1/3 + the three batch-19 fixes; lstol 0: t = 1
with 843 failed attempts):
  L8D_ls      + line search (lstol 0 -> 0.9): the like-for-like control
  L8D_fitall  L8D_ls + every non-P-arcus loft fitted: task 2's combined model
  L8D_conn    L8D_ls + all 12 non-P-arcus lofts replaced by their Abaqus connectors as springs: the faithful reference
              for the Abaqus comparison
usage: py -3.10 build_batch22.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit, LOFT_FAMILIES, NEAR_ZERO_LOFTS
from loft_fit_all import fit_all

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
LS = os.path.join(RUNS, 'L6_ls', 'L6_ls.feb')
DX = os.path.join(RUNS, 'L6D_x', 'L6D_x.feb')
TISSUE = 1.06e-9
LOFTS = ('AVW-Para-L_fan', 'AVW-Para-R_fan', 'CL-L_fan', 'CL-R_fan', 'USL-L_fan', 'USL-R_fan', 'PM_fan') + NEAR_ZERO_LOFTS
WANT = set(sys.argv[1:])
FITS = {}


def want(name):
    return not WANT or name in WANT


def new(name, base):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(base)


def fit(m, base, lofts):
    if base not in FITS:
        FITS[base] = fit_all(base, list(LOFTS))
    for loft in lofts:
        r = FITS[base][loft]
        m.set_loft_material(loft, r['c'], r['m'], 250.0, density=TISSUE if loft in NEAR_ZERO_LOFTS else None,
                            note=f"loft_fit_all.py strip fit to its {r['n']} Abaqus connectors, rms {100 * r['err']:.1f} % "
                                 f"over u = 2-47 mm")


if want('L8_fitall'):
    m = new('L8_fitall', LS); fit(m, LS, LOFTS); emit('L8_fitall', m)
if want('L8_nzweak'):
    m = new('L8_nzweak', LS); fit(m, LS, NEAR_ZERO_LOFTS); emit('L8_nzweak', m)
if want('L8_aggr'):
    m = new('L8_aggr', LS); m.set_solver(aggressiveness=1, cutback=0.5); emit('L8_aggr', m)
if want('L8D_ls'):
    m = new('L8D_ls', DX); m.set_solver(lstol=0.9); emit('L8D_ls', m)
if want('L8D_fitall'):
    m = new('L8D_fitall', DX); m.set_solver(lstol=0.9); fit(m, DX, LOFTS); emit('L8D_fitall', m)
if want('L8D_conn'):
    m = new('L8D_conn', DX); m.set_solver(lstol=0.9); m.lofts_to_connectors(LOFT_FAMILIES); emit('L8D_conn', m)
