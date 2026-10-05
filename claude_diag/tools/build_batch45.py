"""Batch 45 (2026-09-26): connective-tissue changes for the perineal body's rotation about -x (README_2026-09-26.md).

Batch 44 early readings (t 0.16-0.38): cutting all 26 P-arcus connectors triples the early rotation (springs +6.1 deg at
t 0.17, base +1.8); cutting the posterior sphincter connectors turns the body the wrong way from the start (they hold its
posterior edge up while the PVW takes the anterior edge down); the PVW-LA contact never touches (control identical).
So a real weakening of P-arcus, one change on each L26 line (NOT IN SOURCE: connective tissue, not the vaginal wall):
  L37_{line}_parcus50   Parcus_conn material scale 1 -> 0.5 (half the force at every elongation)
usage: py -3.10 build_batch45.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from build_batch44 import BASES, spring_scale

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

VARIANTS = {'parcus50': lambda m: spring_scale(m, 'Parcus_conn', 0.5, label='NOT IN SOURCE (connective tissue)')}
for line, base in BASES.items():
    for tag, fn in VARIANTS.items():
        name = f'L37_{line}_{tag}'
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(base)
        fn(m)
        emit(name, m)
