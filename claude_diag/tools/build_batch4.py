import os
from variants import Model, JOBS, emit

TRUSS = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.feb')
V32 = os.path.join(JOBS, 'PVP3DModel_v32.feb')
CORR_CTRL = dict(max_retries=20, dtmin=1e-09, cutback=0.25, max_refs=25)
FANS = ['Sphincter-L_fan', 'Sphincter-R_fan', 'Sphincter-P_fan']

# T3f: most Abaqus-faithful LA model: exact ATLA/Arcus ties + 40 sphincter connectors,
# redundant sphincter fans removed, remesh-lost pins restored
m = Model(TRUSS)
m.remove_truss_ties()
m.unmerge_la()
m.add_abaqus_ties_lc()
m.add_all_sphincter_connectors(which=('L', 'R', 'P'))
m.restore_remesh_lost_pins()
m.remove_domains(FANS)
emit('T3f_faithful', m)

# N1f: no-LA reference with the same fan removal + pins
m = Model(V32)
m.set_control(**CORR_CTRL)
m.restore_remesh_lost_pins()
m.remove_domains(FANS)
emit('N1f_v32_pins_nofans', m)
