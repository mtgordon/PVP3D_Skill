"""Batch 15d: LP2DM_m10 + settle phase (dynamic relaxation to the full-load equilibrium)."""
import os
from variants3 import Model3, JOBS, emit

m = Model3(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_m10', 'LP2DM_m10.feb'))
m.add_settle_phase(t_end=2.0, C=20.0)
emit('LP2DM_m10_settle', m)
