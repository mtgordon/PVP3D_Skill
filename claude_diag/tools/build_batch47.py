"""Batch 47 (2026-09-26): the lofts line crawls at t ~0.45 once P-arcus is weakened (L37_lofts_parcus50: 4 failed, the
fastest nodes where the USL-L loft is tied to the PVW, x -16 y -25 z +8, and in the AVW-Para lofts). The USL-L loft is
the one loft still tied by tied-node-on-facet (USL-R shares nodes); the ties take exactly the tissue nodes that sit on loft
nodes, so putting it on shared nodes is equivalent (merge_tnof, as L33_pvw50_clusl50_usltie). One change:
  L39_lofts_parcus50_usltie   L37_lofts_parcus50 + the USL-L loft on shared nodes (a loft convention, NOT IN SOURCE)
usage: py -3.10 build_batch47.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from tnof_merge import merge_tnof, USL_L_TIES

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L39_lofts_parcus50_usltie', 'L37_lofts_parcus50'),)
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        merge_tnof(m, USL_L_TIES)
        emit(name, m)
