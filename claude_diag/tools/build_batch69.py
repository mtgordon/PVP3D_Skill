"""Batch 69 (2026-09-26 ~12:37): the last steps toward 90 deg.

L59_springs_pc0_cu30_antside0 (P-arcus cut, CL/USL 30 %, the 3 anterior pairs of side sphincters cut) peaks +76.8 deg at
t 0.37 with 2 failed attempts; the lofts twin of the recipe + the anterior cut (L59_lofts...) +72.7 at t 0.82. One change
each, NOT IN SOURCE (connective tissue):
  L62_lofts_pc0_cu30_avwparaconn_antside0  L59_lofts_pc10_cu30_avwparaconn_antside0 with P-arcus 10 % -> cut (the lofts
                                           twin of L59_springs_pc0_cu30_antside0)
  L62_springs_pc0_cu30_antside4            L59_springs_pc0_cu30_antside0 + the 4th side pair (body end y 4.7, moment
                                           +0.4 N mm per N: neutral) cut too; the 4 posterior pairs stay
usage: py -3.10 build_batch69.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

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


def side_pair_y47(m):
    X = m.nodes()
    peb = set()
    for blk in m.mesh.findall('Elements'):
        if blk.get('name') == '_PickedSet66':
            for e in blk:
                peb.update(int(v) for v in e.text.split(','))
    ds = next(d for d in m.mesh.findall('DiscreteSet') if d.get('name') == 'LA_sphincter_side_conn')
    cut = []
    for d in list(ds.findall('delem')):
        a, b = (int(v) for v in d.text.split(','))
        t = a if a in peb else b
        if 4.0 <= X[t][1] < 5.5:
            ds.remove(d)
            cut.append(round(float(X[t][1]), 1))
    assert len(cut) == 2, cut
    m.log.append(f'NOT IN SOURCE (connective tissue): the side sphincter pair at body end y {cut[0]} cut as well; '
                 f'{len(ds.findall("delem"))} posterior side connectors kept')


BUILDS = (('L62_lofts_pc0_cu30_avwparaconn_antside0', 'L59_lofts_pc10_cu30_avwparaconn_antside0', parcus_10_to_0),
          ('L62_springs_pc0_cu30_antside4', 'L59_springs_pc0_cu30_antside0', side_pair_y47))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
