"""Batch 20 (2026-09-24): on the combined batch-19 model L5X_all (Load-LA x 1/3 NOT IN SOURCE, tnof ties 100 N/mm,
P-arcus loft edge nodes on the chain, 328 midline stab springs removed), ending at t = 1.0 (full load, no settle),
one change each:
  L6_ctl     the control
  L6_ls      line search on (lstol 0 -> 0.9)            } from the user's '3 Surfaces Traction Working' test
  L6_fn      full Newton (max_ups 10 -> 0)               } (mini3s/: without the line search it failed at 0.46,
  L6_blk     that test's solver block (lstol 0.9, max_ups 0, opt_iter 45, max_refs 100, cutback 0.5)
  L6_pmpeb   the PM_PeB lofts get their own material fitted to their 8 connectors per side (near zero, c1 ~4e-5),
             with tissue density (their connectors are massless; at 7.8e-07 each loft weighed 64 g); and, the user's
             choice at the batch start, the PM_avw_bottom lofts too (2 connectors per side, the same near-zero table;
             c1 ~1.4e-4; ~8 g each at 7.8e-07)
  L6_pebc    (added at the batch start, user) PeB-constrin its own fitted material (26 connectors, the same near-zero
             table; c1 ~6.8e-4), tissue density (at 7.8e-07 the loft weighed ~150 g; it is compressed ~4.6 mm in the
             runs, which the source connectors cannot carry)
  L6_fitall  every non-P-arcus loft its own fitted material (loft_fit_all.py); the near-zero ones (PM_PeB,
             PM_avw_bottom, PeB-constrin) with tissue density
  L6_tie1    (added at the batch start, user; queued behind the first 8) the tied-node-on-facet penalty of the two
             tied lofts (AVW-Para-L, USL-L) 100 -> 1 N/mm, ~10 % of the loft's E t (21 x 0.49 = 10.3 N/mm): the
             3-surface test's lesson (a tie ~1 % of the shell's E t worked there, one as stiff as the shell failed)
Like-for-like base (LPFmsTk_DM1, Abaqus chain mass, to t = 1):
  L6D_x      Load-LA x 1/3 + the three batch-19 fixes (ties, arcus edge, midline springs)
usage: py -3.10 build_batch20.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants5 import Model5, JOBS, emit, remove_stab_springs, NEAR_ZERO_LOFTS
from loft_fit_all import fit_all

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L5X_all', 'L5X_all.feb')
DM1 = os.path.join(RUNS, 'LPFmsTk_DM1', 'LPFmsTk_DM1.feb')
TISSUE = 1.06e-9
LOFTS = ('AVW-Para-L_fan', 'AVW-Para-R_fan', 'CL-L_fan', 'CL-R_fan', 'USL-L_fan', 'USL-R_fan', 'PM_fan') + NEAR_ZERO_LOFTS
PM_PEB = ('PM_PeB_Left_fan', 'PM_PeB_Right_fan')
PM_AVW = ('PM_avw_bottom_left_fan', 'PM_avw_bottom_right_fan')
WANT = set(sys.argv[1:])

fits = fit_all(BASE, list(LOFTS))


def want(name):
    return not WANT or name in WANT


def new(name, base=BASE, t_end=1.0):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model5(base)
    if t_end:
        m.set_end_time(t_end)
    return m


def fit(m, lofts):
    for loft in lofts:
        r = fits[loft]
        m.set_loft_material(loft, r['c'], r['m'], 250.0, density=TISSUE if loft in NEAR_ZERO_LOFTS else None,
                            note=f"loft_fit_all.py strip fit to its {r['n']} Abaqus connectors, rms {100 * r['err']:.1f} % "
                                 f"over u = 2-47 mm")


if want('L6_ctl'):
    emit('L6_ctl', new('L6_ctl'))
if want('L6_ls'):
    m = new('L6_ls'); m.set_solver(lstol=0.9); emit('L6_ls', m)
if want('L6_fn'):
    m = new('L6_fn'); m.set_solver(max_ups=0); emit('L6_fn', m)
if want('L6_blk'):
    m = new('L6_blk'); m.set_solver(lstol=0.9, max_ups=0, opt_iter=45, max_refs=100, cutback=0.5); emit('L6_blk', m)
if want('L6_pmpeb'):
    m = new('L6_pmpeb'); fit(m, PM_PEB + PM_AVW); emit('L6_pmpeb', m)
if want('L6_pebc'):
    m = new('L6_pebc'); fit(m, ('PeB-constrin_fan',)); emit('L6_pebc', m)
if want('L6_fitall'):
    m = new('L6_fitall'); fit(m, LOFTS); emit('L6_fitall', m)
if want('L6_tie1'):
    m = new('L6_tie1'); m.set_tnof_penalty(1.0); emit('L6_tie1', m)
if want('L6D_x'):
    m = new('L6D_x', base=DM1, t_end=None)
    m.scale_surface_load('Load-LA', 1 / 3)
    m.set_tnof_penalty(100.0)
    m.parcus_edge_to_chain()
    remove_stab_springs(m, only_nodeset='BC-VW-mid')
    emit('L6D_x', m)
