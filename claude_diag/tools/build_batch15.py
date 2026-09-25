"""Batch 15: Load-LA on TRUSS_fixed2, then LA material fidelity, then two Load-LA strategies (one change per row).

LP0  TRUSS_fixed2 + Abaqus Load-LA                     (baseline for this batch)
LP1  LP0 + LA Yeoh (I1-only, Marlow-equivalent), k from mu0
LP2  LP1 + LA k = 1.0 (Abaqus constant-nu volume response over 20-50 % strain)
LP2D LP2 + implicit DYNAMIC, Abaqus densities (source is *Dynamic, Explicit)
LP2R LP2 + Load-LA ramped on its own curve from t = 0.5 to 1.0
"""
from variants3 import Model3, TRUSS_FIXED2, emit


def lp(*ops):
    m = Model3(TRUSS_FIXED2)
    m.add_la_pressure()
    for op in ops:
        op(m)
    return m


yeoh = lambda m: m.set_la_yeoh()
k1 = lambda m: m.set_la_k(1.0)
emit('LP0_laload', lp())
emit('LP1_yeoh', lp(yeoh))
emit('LP2_yeoh_k1', lp(yeoh, k1))
emit('LP2D_yeoh_k1_dyn', lp(yeoh, k1, lambda m: m.set_dynamic(rhoi=0.5)))
emit('LP2R_yeoh_k1_lalate', lp(yeoh, k1, lambda m: m.delay_la_load(0.5, 1.0)))
