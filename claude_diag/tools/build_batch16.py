"""Batch 16: the Abaqus posterior-arcus -> vaginal-wall link (P-arcus-L/R-1..13), one change per run.
Static on LP2 (the chain-buckling wall at t ~ 0.29) and dynamic + settle on LP2DM_m10_settle (full-load equilibrium)."""
import os
from variants4 import Model4, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
LP2 = os.path.join(RUNS, 'LP2_yeoh_k1', 'LP2_yeoh_k1.feb')
SETTLE = os.path.join(RUNS, 'LP2DM_m10_settle', 'LP2DM_m10_settle.feb')

for name, base, op in (('LPF_fanties', LP2, 'tie_parcus_fans'),
                       ('LPC_pconn', LP2, 'add_parcus_connectors'),
                       ('LPF_DMs_fanties', SETTLE, 'tie_parcus_fans'),
                       ('LPC_DMs_pconn', SETTLE, 'add_parcus_connectors')):
    m = Model4(base)
    getattr(m, op)()
    emit(name, m)
