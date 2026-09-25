"""Batch 16g. (1) LPBp_beam_pinned with a non-symmetric global stiffness (a geometrically exact beam's tangent is
not symmetric away from equilibrium; the model runs with symmetric_stiffness = preferred). (2) The P-arcus fans tied
on the faithful-mass dynamic run, like-for-like with Abaqus/Explicit (pairs with LPC_DM1_pconn)."""
import os
from variants4 import Model4, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
m = Model4(os.path.join(RUNS, 'LPBp_beam_pinned', 'LPBp_beam_pinned.feb'))
s = m.root.find('Control').find('solver').find('symmetric_stiffness')
old = s.text
s.text = 'non-symmetric'
m.log.append(f'solver symmetric_stiffness {old} -> non-symmetric')
emit('LPBpN_beam_pinned_nonsym', m)

m = Model4(os.path.join(RUNS, 'LP2DM_dyn_chainmass', 'LP2DM_dyn_chainmass.feb'))
m.tie_parcus_fans()
emit('LPF_DM1_fanties', m)
