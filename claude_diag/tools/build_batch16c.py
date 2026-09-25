"""Batch 16c: diagnostic rerun of LPF_fanties (fans tied, static) to localize its t = 0.5705 failure:
per-iteration plot (displacement + relative volume only), stop at t = 0.58, run single-threaded with the
console's negative-jacobian listing (see run command in README_2026-09-23c.md)."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LPF_fanties', 'LPF_fanties.feb'))
m.set_diagnostic_output()
m.set_control(time_steps=58)
emit('LPF_diag', m)
