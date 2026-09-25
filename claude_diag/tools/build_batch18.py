"""Batch 18 (2026-09-24): every loft against its own Abaqus connectors (README_2026-09-24.md), one change per variant,
on the fast settle base LPFmsTk_m10s (chain mass x10, merged + fitted P-arcus lofts, settle damping from t = 1),
screened to t = 1.2 (run length only). L4_ctl is the base itself to t = 1.2, run alongside for a fair comparison
of retries and wall time."""
import os

from variants5 import Model5, JOBS, emit
from loft_fit_all import fit_all

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'LPFmsTk_m10s', 'LPFmsTk_m10s.feb')
T_END = 1.2
FITTED = {'L4A_avwpara': ('AVW-Para-L_fan', 'AVW-Para-R_fan'), 'L4C_cl': ('CL-L_fan', 'CL-R_fan'),
          'L4U_usl': ('USL-L_fan', 'USL-R_fan'), 'L4P_pm': ('PM_fan',)}

fits = fit_all(BASE, [l for ls in FITTED.values() for l in ls])


def new(name):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model5(BASE)
    m.set_end_time(T_END)
    return m


emit('L4_ctl', new('L4_ctl'))
for name, lofts in FITTED.items():
    m = new(name)
    for loft in lofts:
        r = fits[loft]
        m.set_loft_material(loft, r['c'], r['m'], 250.0,
                            note=f"loft_fit_all.py strip fit to its {r['n']} Abaqus connectors, rms {100 * r['err']:.1f} % "
                                 f"over u = 2-47 mm")
    emit(name, m)
m = new('L4R_rm5')
m.remove_lofts()
emit('L4R_rm5', m)
m = new('L4D_rho')
m.set_loft_density(1.06e-9)
emit('L4D_rho', m)
