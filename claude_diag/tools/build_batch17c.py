"""Batch 17c: the loft configuration (merged onto the arcus nodes, soft Ogden fitted to the connectors, 0.49 mm) on the
fast base LP2DM_m10_settle (chain mass x10, loads held t = 1 -> 2 with mass damping C = 20/s): the full-load
equilibrium quickly (the masses change the path, not the equilibrium). Run with 8 threads and -dump=1."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_m10_settle', 'LP2DM_m10_settle.feb'))
m.truss2_to_line2(); m.merge_parcus_fans(); m.set_parcus_fan_material(); m.set_parcus_fan_thickness(0.49)
emit('LPFmsT_m10s', m)
