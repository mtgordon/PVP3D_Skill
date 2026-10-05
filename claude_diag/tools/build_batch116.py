"""Batch 116 (2026-09-29): the all-lofts edge holds after the bug fix, as penalty springs (penalty 100 N/mm, maxaug 0)
instead of augmented ties. L115_fibre_all_m05_edges (augmented, fixed) steps at dt ~2e-4 with 2-4 augmentations a step
(the P-arcus holds outside tolerance). One change from L115_fibre_all_m05_edges. NOT IN SOURCE.
  L116_fibre_all_m05_edgespen
usage: py -3.10 build_batch116.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_batch109 as b9  # noqa: E402
from build_batch109 import Model109, fibre, ALL, BASE, RUNS, emit  # noqa: E402

b9.HOLD.update(penalty='100', maxaug='0')
name = 'L116_fibre_all_m05_edgespen'
if __name__ == '__main__':
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    mdl = Model109(os.path.join(RUNS, BASE, BASE + '.feb'))
    fibre(ALL, 0.05, edges=True)(mdl)
    emit(name, mdl)
