"""Batch 17d: the soft Ogden loft swells (J median 1.32, max 4.25 by t = 0.96) because k = 50 mu0 is small once
the m1 = 3.5 law has stiffened ~5x at ~3x stretch. One change from LPFmsT_DM1: loft k = 250 mu0 (nu ~ 0.49 at the
working stretch; the loft then thins like the incompressible strip fit assumed)."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_dyn_chainmass', 'LP2DM_dyn_chainmass.feb'))
m.truss2_to_line2(); m.merge_parcus_fans(); m.set_parcus_fan_material(k_over_mu=250.0); m.set_parcus_fan_thickness(0.49)
emit('LPFmsTk_DM1', m)
