"""Batch 16d: NOT IN SOURCE numerical stabiliser - a linear-beam along the posterior-arcus springs (bending/shear
stiffness, ~0 axial), one change on each static base: LP2 (no arcus link), LPF_fanties (fans tied), LPC_pconn
(26 connectors)."""
import os
from variants4 import Model4, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
for name, base in (('LPB_beam', 'LP2_yeoh_k1'), ('LPFB_fanties_beam', 'LPF_fanties'), ('LPCB_pconn_beam', 'LPC_pconn')):
    m = Model4(os.path.join(RUNS, base, base + '.feb'))
    m.add_arcus_beam()
    emit(name, m)
