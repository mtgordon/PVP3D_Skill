"""Batch 83 (2026-09-27 ~17:30): the Yeoh walls on the lines that now run crawl at two contact edges.

With PVW_LA seg_up 2 the faithful (Yeoh) vaginal walls get further (L78_nopinch_pvwsegup2_vwyeoh: t 0.57 in 50 min at the
source loads; without it L74_nopinch_vwyeoh took 7.5 h to t 0.46), but on the springs line at Load-LA 1/3
(L78_springs_la3_vwyeoh, t 0.40) two things rattle: the PVW's lateral edge row on the LA (1983/1984, 1899, 5402, 7327 at
2-4 m/s, 0.01-0.03 MPa; the softer walls bring the PVW edge onto the LA with some pressure) and the canal seam (PVW edge
nodes 2574-2576, 2347 at ~0.7 m/s). One change each:
  L79_springs_la3_vwyeoh_pvwpen05    L78_springs_la3_vwyeoh, PVW_LA penalty 5 -> 0.5 (seg_up 2 and the 30 facets kept)
  L79_springs_la3_vwyeoh_se1pen05    L78_springs_la3_vwyeoh, SlidingElastic1 penalty 5 -> 0.5 (the canal)
  L79_springs_la3_vwyeoh_pvwedge6    L78_springs_la3_vwyeoh, NOT IN SOURCE: the 26 PVW_LA primary facets holding an edge
                                     node within 6 mm of an LA node at rest removed (covers the flipping edge nodes, all
                                     4.6-5.4 mm from the LA; 620 of 646 left)
  L79_nopinch_pvwsegup2_vwyeoh_se1pen05  L78_nopinch_pvwsegup2_vwyeoh (source loads), SlidingElastic1 penalty 5 -> 0.5
usage: py -3.10 build_batch83.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import collections
import os
import sys

import numpy as np

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def contact(name, **kw):
    def fn(m):
        assert m.set_contact(name, **kw) == 1
    return fn


def pvwedge(m, tol=6.0):
    X = m.nodes()
    surf = {s.get('name'): s for s in m.mesh.findall('Surface')}
    P = surf['PVW_LA_primary']
    fl = lambda s: [tuple(int(v) for v in f.text.split(',')) for f in s]
    c = collections.Counter()
    for f in fl(P):
        for i in range(len(f)):
            c[tuple(sorted((f[i], f[(i + 1) % len(f)])))] += 1
    bd = {v for e, k in c.items() if k == 1 for v in e}
    LA = np.array([X[v] for v in {v for f in fl(surf['PVW_LA_secondary']) for v in f}])
    near = {v for v in bd if np.min(np.linalg.norm(LA - X[v], axis=1)) < tol}
    n0 = len(P)
    gone = [f for f in list(P) if {int(v) for v in f.text.split(',')} & near]
    for f in gone:
        P.remove(f)
    for i, f in enumerate(P, 1):
        f.set('id', str(i))
    m.log.append(f'NOT IN SOURCE (diagnostic): {len(gone)} PVW_LA primary facets removed: those holding one of the '
                 f'{len(near)} edge nodes within {tol} mm of an LA node at rest (the PVW edge rows that meet the LA); '
                 f'{len(P)} of {n0} left')


BUILDS = (('L79_springs_la3_vwyeoh_pvwpen05', 'L78_springs_la3_vwyeoh', contact('PVW_LA', penalty=0.5)),
          ('L79_springs_la3_vwyeoh_se1pen05', 'L78_springs_la3_vwyeoh', contact('SlidingElastic1', penalty=0.5)),
          ('L79_springs_la3_vwyeoh_pvwedge6', 'L78_springs_la3_vwyeoh', pvwedge),
          ('L79_nopinch_pvwsegup2_vwyeoh_se1pen05', 'L78_nopinch_pvwsegup2_vwyeoh',
           contact('SlidingElastic1', penalty=0.5)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
