import os
from variants import Model, JOBS
from variants import emit

TRUSS = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.feb')
BASE = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_CORRECTED.feb')

# T1: TRUSS with the phantom-ribbon ties replaced by the Abaqus ties as exact linear constraints
m = Model(TRUSS)
m.remove_truss_ties()
m.add_abaqus_ties_lc()
emit('T1_truss_lc', m)

# T2: T1 + LA unmerged from the L/R fans, replaced by the 16 Abaqus L/R sphincter connectors
m = Model(TRUSS)
m.remove_truss_ties()
m.unmerge_la()
m.add_abaqus_ties_lc()
m.add_all_sphincter_connectors(which=('L', 'R'))
emit('T2_truss_lc_connLR', m)

# T3: T2 + the 24 posterior sphincter connectors (all 40 Abaqus LA-PeB connectors)
m = Model(TRUSS)
m.remove_truss_ties()
m.unmerge_la()
m.add_abaqus_ties_lc()
m.add_all_sphincter_connectors(which=('L', 'R', 'P'))
emit('T3_truss_lc_connLRP', m)

# E5: tube model (no LA-truss ties, as in your 0.219 run) + all 40 connectors instead of the fan merge
m = Model(BASE)
m.unmerge_la()
m.add_all_sphincter_connectors(which=('L', 'R', 'P'))
emit('E5_tubes_connLRP', m)
