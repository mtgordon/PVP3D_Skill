"""Batch 14: stack the Abaqus-faithful fixes that each ran on B1 (combination check)."""
import os
from variants2 import Model2, JOBS, emit, set_shell_normal_nodal, FANS
from build_batch10 import peb_bottom_facets

B1P1 = os.path.join(JOBS, 'claude_diag', 'runs', 'B1P1_vwmid_pvwla', 'B1P1_vwmid_pvwla.feb')
fac = peb_bottom_facets()
# ALL1 = B1 + P1 (contact) + C1 (PeB facets) + X1 (connectors constant)
m = Model2(B1P1); m.add_peb_bottom_to_avw_contact(fac); m.set_connector_extend('constant'); emit('ALL1_faithful', m)
# ALL2 = ALL1 + F1 (all fans element normals)
m = Model2(B1P1); m.add_peb_bottom_to_avw_contact(fac); m.set_connector_extend('constant'); set_shell_normal_nodal(m, FANS)
emit('ALL2_faithful_fans', m)
