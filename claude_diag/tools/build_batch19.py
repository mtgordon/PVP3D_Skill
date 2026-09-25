"""Batch 19 (2026-09-24, after the user's review): the LA load at 1/3 (NOT IN SOURCE; the user's choice while the LA
deforms far more than expected of Abaqus), plus the attachment and midline fixes, one change each.

Fast base (LPFmsTk_m10s: chain mass x10 + settle, to t = 1.2):
  L5_ctl    Load-LA x 1/3 (the new control)
  L5T_ties  + the two tied-node-on-facet lofts (AVW-Para-L, USL-L) tied properly: penalty 0.0005 -> 100 N/mm
  L5G_gaps  + the 8 loose P-arcus loft edge nodes per side made arcus chain nodes (chain springs split)
  L5S_stab  + the 328 midline zero-length stab_* ground springs removed (not in the source)
  L5X_all   + all three
Like-for-like base (LPFmsTk_DM1: Abaqus chain mass, 1 s ramp):
  L5D_la3   Load-LA x 1/3 only
"""
import os

from variants5 import Model5, JOBS, emit, remove_stab_springs

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
FAST = os.path.join(RUNS, 'LPFmsTk_m10s', 'LPFmsTk_m10s.feb')
DM1 = os.path.join(RUNS, 'LPFmsTk_DM1', 'LPFmsTk_DM1.feb')


def new(name, base=FAST, t_end=1.2):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model5(base)
    if t_end:
        m.set_end_time(t_end)
    m.scale_surface_load('Load-LA', 1 / 3)
    return m


def ties(m):
    m.set_tnof_penalty(100.0)


def gaps(m):
    m.parcus_edge_to_chain()


def stab(m):
    remove_stab_springs(m, only_nodeset='BC-VW-mid')


emit('L5_ctl', new('L5_ctl'))
for name, ops in (('L5T_ties', (ties,)), ('L5G_gaps', (gaps,)), ('L5S_stab', (stab,)), ('L5X_all', (ties, gaps, stab))):
    m = new(name)
    for op in ops:
        op(m)
    emit(name, m)
emit('L5D_la3', new('L5D_la3', base=DM1, t_end=None))
