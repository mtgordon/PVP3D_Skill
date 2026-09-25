"""Batch 15c (one change per row).

LP2S3   LP2 + LA x3 stiffer than Abaqus (static; user-approved non-faithful test)
LP2S10  LP2 + LA x10 stiffer than Abaqus (static)
LP2DM_r0   LP2DM + rhoi 0 (maximum high-frequency damping instead of 0.5)
LP2DM_m10  LP2DM with the chain mass x10 (mass scaling on the chains only)
"""
from variants3 import Model3, TRUSS_FIXED2, emit


def lp2(*ops):
    m = Model3(TRUSS_FIXED2)
    m.add_la_pressure()
    m.set_la_yeoh()
    m.set_la_k(1.0)
    for op in ops:
        op(m)
    return m


emit('LP2S3_la_x3', lp2(lambda m: m.scale_la(3)))
emit('LP2S10_la_x10', lp2(lambda m: m.scale_la(10)))
emit('LP2DM_r0', lp2(lambda m: m.add_chain_mass(), lambda m: m.set_dynamic(rhoi=0)))
emit('LP2DM_m10', lp2(lambda m: m.add_chain_mass(), lambda m: m.set_dynamic(rhoi=0.5),
                      lambda m: m.log.append('chain mass x10 (mass scaling on the chains only)'),
                      lambda m: [e.find('density').__setattr__('text', '0.0011')
                                 for e in m.root.find('Material') if e.get('name') == 'chain_mass_mat']))
