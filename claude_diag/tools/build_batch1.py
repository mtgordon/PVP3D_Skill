import os
import shutil
from variants import Model, JOBS

BASE = os.path.join(JOBS, 'PVP3DModel_v21_v3_HYBRID_withLA_CORRECTED.feb')
V32 = os.path.join(JOBS, 'PVP3DModel_v32.feb')
OUT = os.path.join(JOBS, 'claude_diag', 'runs')

CORR_CTRL = dict(max_retries=20, dtmin=1e-09, cutback=0.25, max_refs=25)
V32_CTRL = dict(max_retries=8, dtmin=1e-06, cutback=0.5, max_refs=12)


def emit(name, model=None, src=None):
    d = os.path.join(OUT, name)
    os.makedirs(d, exist_ok=True)
    dst = os.path.join(d, name + '.feb')
    if model is None:
        shutil.copy(src, dst)
    else:
        model.write(dst)
    print('wrote', dst)


# E0: control, byte-identical copy of the base
emit('E0_base', src=BASE)

# E1: LA unmerged from the Sphincter-L/R fans (LA held only by its origin BCs)
m = Model(BASE)
m.unmerge_la()
emit('E1_unmerge', m)

# E2: E1 + the 8 Abaqus LA->PeB sphincter connectors as nonlinear springs
m = Model(BASE)
m.unmerge_la()
m.add_sphincter_connectors()
emit('E2_abqconn', m)

# E3: no-LA v32 with the CORRECTED solver settings (isolates solver-setting effect)
m = Model(V32)
m.set_control(**CORR_CTRL)
emit('E3_v32_corrctrl', m)

# E4: base with the v32 solver settings
m = Model(BASE)
m.set_control(**V32_CTRL)
emit('E4_base_v32ctrl', m)
