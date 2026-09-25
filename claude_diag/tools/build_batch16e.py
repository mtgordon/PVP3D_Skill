"""Batch 16e: diagnostic rerun of LPB_beam (early failure at t ~ 0.022): per-iteration plot, stop at t = 0.025."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LPB_beam', 'LPB_beam.feb'))
m.set_diagnostic_output()
m.set_control(time_steps=5)
emit('LPB_diag', m)
