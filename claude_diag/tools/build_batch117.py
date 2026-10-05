"""Batch 117 (2026-09-29): the all-lofts edge-hold fibre lofts with the stiffening fibre capped. L115_fibre_all_m05_edges
and L116_fibre_all_m05_edgespen both died at t ~0.545: single USL / AVW-Para loft elements stretched 2-3.8 along their
fibres, beyond the laws' fitted range, where fiber-exp-pow's exponential / power term gave 73-3640 MPa. Here the
stiffening fibre is fiber-exp-pow-linear (FEBioMech/FEFiberPowLinear.cpp; checked on one element against its formula):
exactly the fitted exp-pow below lam0 = 1 + u_max / (the family's shortest connector), linear above (stress and slope
continuous). One change from L116_fibre_all_m05_edgespen (held edges as penalty springs 100 N/mm). NOT IN SOURCE.
  L117_fibre_all_m05_edges_cap
usage: py -3.10 build_batch117.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_batch109 as b9  # noqa: E402
from build_batch109 import Model109, fibre, ALL, BASE, RUNS, emit  # noqa: E402

b9.HOLD.update(penalty='100', maxaug='0')
b9.CAP['on'] = True
name = 'L117_fibre_all_m05_edges_cap'
if __name__ == '__main__':
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    mdl = Model109(os.path.join(RUNS, BASE, BASE + '.feb'))
    fibre(ALL, 0.05, edges=True)(mdl)
    emit(name, mdl)
