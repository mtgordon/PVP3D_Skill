"""Batch 12: on top of B1 (Abaqus BC-VW-mid XSYMM restored)."""
import os
from variants2 import Model2, JOBS, emit, remove_stab_springs

B1 = os.path.join(JOBS, 'claude_diag', 'runs', 'B1_vwmid', 'B1_vwmid.feb')
m = Model2(B1); m.add_pvw_la_contact(); emit('B1P1_vwmid_pvwla', m)
m = Model2(B1); remove_stab_springs(m, 'BC-VW-mid'); emit('B2_vwmid_nomidsprings', m)
