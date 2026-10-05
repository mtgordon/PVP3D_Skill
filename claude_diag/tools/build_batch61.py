"""Batch 61 (2026-09-26 ~10:48): the paper's Figure 3 case (springs, L53_springs_pc10_cu30) under the paper's conditions.
One change each from L53_springs_pc10_cu30:
  L54_springs_pc10_cu30_m1      the source arcus chain mass (the paper's Abaqus model has its real masses)
  L54_springs_pc10_cu30_LAfull  Load-LA 0.00467 -> the source 0.014 (the paper applies the full pressure to the levator
                                too); a comparison with the figure's load, not the LA question (the lines stay at 1/3)
usage: py -3.10 build_batch61.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch57 import mass1

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def la_full(m):
    sl = next(s for s in m.root.find('Loads') if s.get('name') == 'Load-LA')
    old = sl.find('pressure').text
    sl.find('pressure').text = '0.014'
    m.log.append(f'Load-LA pressure {old} -> 0.014 MPa (the source value, as the paper\'s Figure 3 case; the lines use 1/3)')


BUILDS = (('L54_springs_pc10_cu30_m1', mass1), ('L54_springs_pc10_cu30_LAfull', la_full))
if __name__ == '__main__':
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, 'L53_springs_pc10_cu30', 'L53_springs_pc10_cu30.feb'))
        fn(m)
        emit(name, m)
