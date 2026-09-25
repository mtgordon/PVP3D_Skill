import os
from variants import Model, JOBS, emit

TRUSS = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.feb')
V32 = os.path.join(JOBS, 'PVP3DModel_v32.feb')
CORR_CTRL = dict(max_retries=20, dtmin=1e-09, cutback=0.25, max_refs=25)

# N1: no-LA v32 (+ same solver settings as the LA runs) with the remesh-lost pins restored
m = Model(V32)
m.set_control(**CORR_CTRL)
m.restore_remesh_lost_pins()
emit('N1_v32_pins', m)

# T2s / T3s: batch-2 T2 / T3 plus the restored pins
for tag, which in (('T2s_truss_lc_connLR_pins', ('L', 'R')), ('T3s_truss_lc_connLRP_pins', ('L', 'R', 'P'))):
    m = Model(TRUSS)
    m.remove_truss_ties()
    m.unmerge_la()
    m.add_abaqus_ties_lc()
    m.add_all_sphincter_connectors(which=which)
    m.restore_remesh_lost_pins()
    emit(tag, m)
