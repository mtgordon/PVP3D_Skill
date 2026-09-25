"""Batch 16b: the 26 Abaqus P-arcus connectors on the faithful-mass dynamic run (like-for-like with Abaqus/Explicit:
1 s smooth-step ramp, Abaqus densities, truss *Density 0.00011 on the chains), for the t = 1 comparison."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_dyn_chainmass', 'LP2DM_dyn_chainmass.feb'))
m.add_parcus_connectors()
emit('LPC_DM1_pconn', m)
