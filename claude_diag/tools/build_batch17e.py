"""Batch 17e: the loft configuration on the fast settle base, with a static finish. One change from LPFmsT_m10s:
the dynamic ramp + settle stops at t = 1.2 (already relaxed there: LA max 52 mm, ~34 mm/s) and a STATIC step follows
with all loads held, to reach the exact full-load static equilibrium."""
import os
from variants4 import Model4, JOBS, emit

m = Model4(os.path.join(JOBS, 'claude_diag', 'runs', 'LP2DM_m10_settle', 'LP2DM_m10_settle.feb'))
m.truss2_to_line2(); m.merge_parcus_fans(); m.set_parcus_fan_material(); m.set_parcus_fan_thickness(0.49)
m.add_static_finish(t_dyn_end=1.2, static_steps=20, static_dt=0.05)
emit('LPFmsT_m10s2', m)
