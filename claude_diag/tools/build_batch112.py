"""Batch 112 (2026-09-29): the all-lofts edge-hold test of batch 110 with the holds as stiff penalty springs.
L110_fibre_all_m05_edges crawled from the start (dt ~2e-4, 2-5 augmentations and repeated stiffness reforms a step;
stopped at t 0.063): the ~200 new edge-hold linear constraints were augmented (maxaug 10) like LA_truss_ties. Here the
same holds with maxaug 0 and penalty 100 N/mm: each a zero-length spring of 100 N/mm per dof between the edge node and
its interpolated point (skill gotcha 27; ~0.01 mm offset at 1 N). A solver setting of the NOT-IN-SOURCE holds, one
change from L110_fibre_all_m05_edges.
  L112_fibre_all_m05_edgespen   L110_fibre_all_m05_edges, edge holds penalty 100 / maxaug 0
usage: py -3.10 build_batch112.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_batch109 as b9  # noqa: E402
from build_batch109 import Model109, fibre, ALL, BASE, RUNS, emit  # noqa: E402

b9.HOLD.update(penalty='100', maxaug='0')
name = 'L112_fibre_all_m05_edgespen'
if __name__ == '__main__':
    if os.path.exists(os.path.join(RUNS, name)):
        print(f'{name} exists; not overwriting')
    else:
        mdl = Model109(os.path.join(RUNS, BASE, BASE + '.feb'))
        fibre(ALL, 0.05, edges=True)(mdl)
        emit(name, mdl)
