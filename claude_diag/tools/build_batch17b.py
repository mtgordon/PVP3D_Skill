"""Batch 17b: LPFms_DM1 crawls at t ~ 0.63 with the soft 0.125 mm loft in a tension field (stretched ~3x along the
connector direction, laterally compressed in ~half its elements -> wrinkling). One change: P-arcus loft thickness
0.125 -> 0.49 mm (as the other lofts) with c1, k and density scaled so the membrane pull and the loft mass are unchanged."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_dyn_chainmass', 'LP2DM_dyn_chainmass.feb'))
m.truss2_to_line2(); m.merge_parcus_fans(); m.set_parcus_fan_material(); m.set_parcus_fan_thickness(0.49)
emit('LPFmsT_DM1', m)
