"""Batch 53 (2026-09-26 ~09:25): the lofts line with the AVW-Para lofts as their connectors.

L44_lofts_parcus25_avwparaconn (L38_lofts_parcus25 + the AVW-Para-L/R lofts replaced by their Abaqus connectors) passes the
t ~0.43 crawl with 1 failed attempt (t 0.62 in 16 min, rot_-x +24.3), where the loft versions crawl (parcus25, the AVW-Para-L
loft on shared nodes, dtol 0.01): the AVW-Para lofts, which can carry compression, are what goes floppy once the PVW's
lateral hold is weakened. Whether the lofts line keeps AVW-Para as a loft is the user's choice; these runs measure it:
  L46_lofts_avwparaconn          L26_lofts_la3_pm + the AVW-Para lofts as connectors (the control: does the swap alone
                                 change the body's rotation?)
  L46_lofts_parcus10_avwparaconn L44_lofts_parcus25_avwparaconn with P-arcus 25 % -> 10 %   NOT IN SOURCE (connective tissue)
usage: py -3.10 build_batch53.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def avwpara_conn(m):
    m.lofts_to_connectors(('AVW-Para-L', 'AVW-Para-R'))
    m.log.append('the AVW-Para-L/R lofts replaced by their connectors (a loft-vs-connector choice for the user)')


def parcus_25_to_10(m):
    disc = m.root.find('Discrete')
    d = next(e for e in disc.findall('discrete') if e.get('discrete_set') == 'Parcus_conn')
    dm = disc.findall('discrete_material')[int(d.get('dmat')) - 1]
    el = dm.find('scale')
    assert abs(float(el.text) - 0.25) < 1e-9, el.text
    el.text = '0.1'
    m.log.append('NOT IN SOURCE (connective tissue): Parcus_conn material scale 0.25 -> 0.1')


BUILDS = (('L46_lofts_avwparaconn', 'L26_lofts_la3_pm', avwpara_conn),
          ('L46_lofts_parcus10_avwparaconn', 'L44_lofts_parcus25_avwparaconn', parcus_25_to_10))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
