"""Batch 63 (2026-09-26 ~11:02): beyond the paper's recipe, toward 90 deg (springs line).

L53_springs_pc10_cu30 (P-arcus 10 %, CL/USL 30 %: the paper's Figure 3 case) turns +44.4 deg by t 0.66. One change each
from it, NOT IN SOURCE (connective tissue):
  L56_springs_pc0_cu30         P-arcus 10 % -> cut (100 % instead of 90 % impaired)
  L56_springs_pc10_cu30_antside0  the 3 anterior pairs of side sphincter connectors cut (they hold the anterior part up;
                               build_batch59.antside0)
usage: py -3.10 build_batch63.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch59 import antside0

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def parcus_10_to_0(m):
    disc = m.root.find('Discrete')
    d = next(e for e in disc.findall('discrete') if e.get('discrete_set') == 'Parcus_conn')
    dm = disc.findall('discrete_material')[int(d.get('dmat')) - 1]
    el = dm.find('scale')
    assert abs(float(el.text) - 0.1) < 1e-9, el.text
    el.text = '0'
    m.log.append('NOT IN SOURCE (connective tissue): Parcus_conn material scale 0.1 -> 0 (all P-arcus connectors cut)')


BUILDS = (('L56_springs_pc0_cu30', parcus_10_to_0), ('L56_springs_pc10_cu30_antside0', antside0))
if __name__ == '__main__':
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, 'L53_springs_pc10_cu30', 'L53_springs_pc10_cu30.feb'))
        fn(m)
        emit(name, m)
