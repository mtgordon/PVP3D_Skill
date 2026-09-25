"""Batch 13: stacking tests, one change each on top of B1 (Abaqus BC-VW-mid XSYMM restored)."""
import os
from variants2 import Model2, JOBS, emit, set_shell_normal_nodal, FANS
from build_batch10 import peb_bottom_facets

B1 = os.path.join(JOBS, 'claude_diag', 'runs', 'B1_vwmid', 'B1_vwmid.feb')
m = Model2(B1); m.add_peb_bottom_to_avw_contact(peb_bottom_facets()); emit('B1C1_vwmid_pebavw', m)
m = Model2(B1); set_shell_normal_nodal(m, FANS); emit('B1F1_vwmid_fans_snn0', m)
m = Model2(B1); m.add_la_pressure(); emit('B1L1_vwmid_laload', m)
m = Model2(B1); m.set_connector_extend('constant'); emit('B1X1_vwmid_extconst', m)
