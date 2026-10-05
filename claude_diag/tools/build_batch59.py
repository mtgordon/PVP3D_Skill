"""Batch 59 (2026-09-26 ~10:35): toward the paper's impaired case (user: panel C, "the part is rotated about 90 degrees and
there is a large PVW bulge"; an impaired case in the paper, so connective-tissue weakening is the route).

Past ~45-55 deg the bulging PVW drags the body forward instead of turning it (L40_springs_noparcus_nosphside: +59.0 at
t 0.65, then dragged 21 mm forward). Moment per N of tension about the body's centroid (+ = the wanted sense), side
sphincter connectors by their body end: y -3.2 .. +2.8 -> +2.9 / +3.3 / +2.7 / +1.7 (they also pull back, -y 0.4-0.6:
against the drag), y 4.7 -> +0.4, y 6.3 / 7.8 / 9.0 -> -0.8 / -1.7 / -2.1 (they hold the anterior part up). So detach only
the 3 anterior pairs (6 springs, body end y >= 5.5), keeping the 10 posterior ones. One change each (a cut = the springs'
<delem> entries removed from LA_sphincter_side_conn), NOT IN SOURCE (connective tissue):
  L52_springs_noparcus_antside0              L36_springs_noparcus (P-arcus cut) + the anterior side sphincters cut
  L52_springs_noparcus_sphside50_antside0    L40_springs_noparcus_sphside50 + the same
  L52_springs_pc10cu50_sphside50_antside0    L48_springs_pc10cu50_sphside50 + the same
  L52_lofts_pc10cu50_sphside50_avwparaconn_antside0  L49_lofts_pc10cu50_sphside50_avwparaconn + the same
usage: py -3.10 build_batch59.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
Y_CUT = 5.5


def antside0(m):
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
        if X[t][1] >= Y_CUT:
            ds.remove(d)
            cut.append(round(float(X[t][1]), 1))
    assert len(cut) == 6, cut
    m.log.append(f'NOT IN SOURCE (connective tissue): the {len(cut)} anterior side sphincter connectors cut (body end y '
                 f'{sorted(set(cut))}; removed from LA_sphincter_side_conn), the {len(ds.findall("delem"))} posterior ones kept')


BUILDS = (('L52_springs_noparcus_antside0', 'L36_springs_noparcus'),
          ('L52_springs_noparcus_sphside50_antside0', 'L40_springs_noparcus_sphside50'),
          ('L52_springs_pc10cu50_sphside50_antside0', 'L48_springs_pc10cu50_sphside50'),
          ('L52_lofts_pc10cu50_sphside50_avwparaconn_antside0', 'L49_lofts_pc10cu50_sphside50_avwparaconn'))
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        antside0(m)
        emit(name, m)
