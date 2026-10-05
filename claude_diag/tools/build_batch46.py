"""Batch 46 (2026-09-26): where and how much the P-arcus connectors hold the perineal body's rotation.

Batch 44 at t ~0.5: cutting all 26 P-arcus connectors (the PVW's lateral edges to the posterior arcus chains) turns the body
+40.7 deg (springs line; base +7.3) and +14.1 at t 0.40 (lofts; base +3.2); the side sphincters are not a restraint (cut:
slightly less rotation); the posterior sphincters are the pivot (cut: the body turns the wrong way). One change each on
both L26 lines:
  L38_{line}_parcus25     Parcus_conn material scale 1 -> 0.25                          NOT IN SOURCE (connective tissue)
  L38_{line}_parcusdist0  DIAGNOSTIC: the 12 distal P-arcus connectors cut (tissue end z <= -22, the lower 6 pairs, from
                          the perineal body's anterior corners up the PVW); the upper 14 kept
  L38_{line}_parcusprox0  DIAGNOSTIC: the 14 proximal ones cut (z >= -18.3); the distal 12 kept
A cut here removes the springs' <delem> entries from the Parcus_conn DiscreteSet (their end nodes stay on the chain and the
PVW elements, so nothing is left orphaned).
usage: py -3.10 build_batch46.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import BASES, spring_scale

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
Z_SPLIT = -20.0      # between the pairs at z -22.3 and -18.3 (tissue ends)


def parcus_cut(m, keep_distal, label='DIAGNOSTIC'):
    X = m.nodes()
    pvw = set()
    for blk in m.mesh.findall('Elements'):
        if blk.get('name') == '_PickedSet64':
            for e in blk:
                pvw.update(int(v) for v in e.text.split(','))
    ds = next(d for d in m.mesh.findall('DiscreteSet') if d.get('name') == 'Parcus_conn')
    cut = []
    for d in list(ds.findall('delem')):
        a, b = (int(v) for v in d.text.split(','))
        t = a if a in pvw else b
        distal = X[t][2] < Z_SPLIT
        if distal != keep_distal:
            ds.remove(d)
            cut.append(round(float(X[t][2]), 1))
    left = len(ds.findall('delem'))
    which = 'proximal' if keep_distal else 'distal'
    m.log.append(f'{label}: the {len(cut)} {which} P-arcus connectors cut (removed from the Parcus_conn DiscreteSet; tissue '
                 f'end z {min(cut)} .. {max(cut)}); {left} kept')


VARIANTS = {'parcus25': lambda m: spring_scale(m, 'Parcus_conn', 0.25, label='NOT IN SOURCE (connective tissue)'),
            'parcusdist0': lambda m: parcus_cut(m, keep_distal=False),
            'parcusprox0': lambda m: parcus_cut(m, keep_distal=True)}
if __name__ == '__main__':
    for line, base in BASES.items():
        for tag, fn in VARIANTS.items():
            name = f'L38_{line}_{tag}'
            if WANT and name not in WANT:
                continue
            assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
            m = Model6(base)
            fn(m)
            emit(name, m)
