"""Batch 51 (2026-09-26 ~09:07): is the lofts line's crawl the AVW-Para lofts carrying compression?

When the proximal P-arcus is weakened the PVW pushes the AVW forward (canal contact); the AVW-Para lofts run from the AVW
to fixed anchors further forward (y ~34), so they may be pushed into compression and wrinkle: a floppy mode a loft can have
and a connector (tension only) cannot. The crawl's fastest nodes are in the AVW and the AVW-Para lofts. One change:
  L44_lofts_parcus25_avwparaconn   L38_lofts_parcus25 + the AVW-Para-L/R lofts replaced by their Abaqus connectors
                                   (lofts_to_connectors, as in the springs line)   DIAGNOSTIC (a loft-vs-connector choice
                                   the user makes)
usage: py -3.10 build_batch51.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
BUILDS = (('L44_lofts_parcus25_avwparaconn', 'L38_lofts_parcus25'),)
if __name__ == '__main__':
    for name, base in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, base, base + '.feb'))
        m.lofts_to_connectors(('AVW-Para-L', 'AVW-Para-R'))
        m.log.append('DIAGNOSTIC: the AVW-Para-L/R lofts replaced by their connectors (is the crawl the lofts in compression?)')
        emit(name, m)
