"""Batch 15e: full-load equilibrium bulge vs LA material (each one change from LP2DM_m10_settle).

LP0DM_m10_settle    LA back to the old 1-term Ogden (c1 0.042712, m1 6.4883, k 0.3488)
LP2S3DM_m10_settle  NOT FAITHFUL: LA Yeoh x3 stiffer (c_i and k x3)
Common: TRUSS_fixed2 + Load-LA + chain truss mass x10 + implicit DYNAMIC (rhoi 0.5) + settle phase to t = 2.
"""
from variants3 import Model3, TRUSS_FIXED2, emit


def dyn_settle(la_ops):
    m = Model3(TRUSS_FIXED2)
    m.add_la_pressure()
    for op in la_ops:
        op(m)
    m.add_chain_mass()
    m.set_dynamic(rhoi=0.5)
    for e in m.root.find('Material'):
        if e.get('name') == 'chain_mass_mat':
            e.find('density').text = '0.0011'
    m.log.append('chain mass x10 (mass scaling on the chains only)')
    m.add_settle_phase(t_end=2.0, C=20.0)
    return m


emit('LP0DM_m10_settle', dyn_settle([]))
emit('LP2S3DM_m10_settle', dyn_settle([lambda m: m.set_la_yeoh(), lambda m: m.set_la_k(1.0), lambda m: m.scale_la(3)]))
