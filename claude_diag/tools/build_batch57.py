"""Batch 57 (2026-09-26 ~10:18): does the fast base's heavier arcus chain bias the P-arcus result?

The P-arcus connectors end on the arcus chain nodes, and the fast base has that chain's mass x10 (chain_mass_mat density
0.0011; the source 0.00011, as L21D/L26D). One change each: the chain density back to the source value.
  L50_springs_la3_pm_m1              L26_springs_la3_pm with the source chain mass
  L50_springs_pc10cu50_m1            L43_springs_parcus10_clusl50 with the source chain mass
  L50_lofts_pc10cu50_avwparaconn_m1  L47_lofts_parcus10_avwparaconn_clusl50 with the source chain mass
usage: py -3.10 build_batch57.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def mass1(m):
    mat = next(x for x in m.root.find('Material') if x.get('name') == 'chain_mass_mat')
    el = mat.find('density')
    old = el.text
    assert abs(float(old) - 0.0011) < 1e-12, old
    el.text = '0.00011'
    m.log.append(f'the arcus chain mass back to the source: chain_mass_mat density {old} -> 0.00011 (the fast base has x10)')


BUILDS = (('L50_springs_la3_pm_m1', 'L26_springs_la3_pm'),
          ('L50_springs_pc10cu50_m1', 'L43_springs_parcus10_clusl50'),
          ('L50_lofts_pc10cu50_avwparaconn_m1', 'L47_lofts_parcus10_avwparaconn_clusl50'))
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        mass1(m)
        emit(name, m)
