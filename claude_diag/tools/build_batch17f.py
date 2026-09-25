"""Batch 17f: one change from LPFmsT_m10s (loft on the fast settle base, which crawled at t ~ 1.258 in the settle
phase with the loft swelling): loft k = 250 mu0 instead of 50 mu0."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_m10_settle', 'LP2DM_m10_settle.feb'))
m.truss2_to_line2(); m.merge_parcus_fans(); m.set_parcus_fan_material(k_over_mu=250.0); m.set_parcus_fan_thickness(0.49)
emit('LPFmsTk_m10s', m)
