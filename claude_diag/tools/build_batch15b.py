"""Batch 15b: LP2D + the Abaqus truss mass on the arcus/ATLA chains (mass-only truss2 elements).

The Load-LA wall is compressive buckling of the massless posterior-arcus spring chains (chain_force.py on
B1F1L1). Abaqus has T3D2 trusses with *Density 0.00011 under *Dynamic, Explicit.
"""
from variants3 import Model3, TRUSS_FIXED2, emit

m = Model3(TRUSS_FIXED2)
m.add_la_pressure()
m.set_la_yeoh()
m.set_la_k(1.0)
m.add_chain_mass()
m.set_dynamic(rhoi=0.5)
emit('LP2DM_dyn_chainmass', m)
