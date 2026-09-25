import os
from variants import Model, JOBS, emit

BASE = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_CORRECTED.feb')
# E0n: the user's CORRECTED model (0.2197) with ONLY the LA element-normal change
m = Model(BASE); m.set_la_shell_normal_nodal(0); emit('E0n_base_snn0', m)
