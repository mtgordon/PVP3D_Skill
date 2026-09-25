import os
from variants import Model, JOBS, emit

TRUSS = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.feb')
FANS = ['Sphincter-L_fan', 'Sphincter-R_fan', 'Sphincter-P_fan']
# Abaqus LA_Yamada50%: Marlow, poisson=0.47. Ogden fit mu0 = c1/2 = 0.021356
K_ABQ = round(2 * 0.021356 * 1.47 / (3 * (1 - 2 * 0.47)), 4)   # 0.3488

def faithful():
    m = Model(TRUSS)
    m.remove_truss_ties()
    m.unmerge_la()
    m.add_abaqus_ties_lc()
    m.add_all_sphincter_connectors(which=('L', 'R', 'P'))
    m.restore_remesh_lost_pins()
    m.remove_domains(FANS)
    return m

if __name__ == "__main__":
  m = faithful(); m.set_la_material(k=K_ABQ); emit('T3f_k', m)
  m = faithful(); m.set_la_shell_type("three-field-shell"); emit("T3f_3f", m)
  m = faithful(); m.set_la_material(k=K_ABQ); m.set_la_shell_type("three-field-shell"); emit("T3f_k3f", m)
  print("K_ABQ", K_ABQ)
