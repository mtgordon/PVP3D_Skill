import os
from variants import Model, JOBS, emit

TRUSS = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.feb')
V32 = os.path.join(JOBS, 'PVP3DModel_v32.feb')
CORR_CTRL = dict(max_retries=20, dtmin=1e-09, cutback=0.25, max_refs=25)
FANS = ['Sphincter-L_fan', 'Sphincter-R_fan', 'Sphincter-P_fan']

def faithful():
    m = Model(TRUSS)
    m.remove_truss_ties()
    m.unmerge_la()
    m.add_abaqus_ties_lc()
    m.add_all_sphincter_connectors(which=('L', 'R', 'P'))
    m.restore_remesh_lost_pins()
    m.remove_domains(FANS)
    return m

# D1: faithful LA model, implicit DYNAMIC with the Abaqus densities (source is *Dynamic, Explicit)
m = faithful()
m.set_dynamic(rhoi=0.5)
emit('D1_faithful_dyn', m)

# D0: no-LA v32 reference, same dynamic settings
m = Model(V32)
m.set_control(**CORR_CTRL)
m.set_dynamic(rhoi=0.5)
emit('D0_v32_dyn', m)
