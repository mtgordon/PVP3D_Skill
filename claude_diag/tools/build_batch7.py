from build_batch6 import faithful, K_ABQ
from variants import emit

m = faithful(); m.set_la_shell_normal_nodal(0); emit('T3f_snn0', m)
m = faithful(); m.set_la_material(k=K_ABQ); m.set_la_shell_normal_nodal(0); emit('T3f_k_snn0', m)
