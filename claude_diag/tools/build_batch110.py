"""Batch 110 (2026-09-29): the fibre lofts with their anchor and tissue edges held along their whole length.
Why: in batch 109 the fibre lofts stretched their connector ends 1.3-3x as far as the springs line's connectors
(fibre_lofts.py --compare). fibre_lofts.loft_strain showed the lofts do not deform along their fibres (P-arcus: fibre
stretch 0.64x its connector line's, the principal stretch > 30 deg off the fibre in 70 % of the elements), and the loft
boundaries are held only at the connector ends: P-arcus 11 of 18 anchor-edge and 13 of 24 tissue-edge nodes, USL-R 13 of
43 tissue-edge nodes, PM_PeB 8 of 18 / 15. Fibres ending at a free edge node carry nothing, which the strip model (every
strip loaded) does not allow for. build_batch109._hold_edges holds every free run of edge nodes between two held nodes of
the same kind (BC anchors: join the BC set; chain or tissue: linear constraints to the two neighbours); the free sides
stay free. NOT IN SOURCE (a loft convention, as variants6.fix_anchor_edge); one change from each batch-109 model.
  L110_fibre_all_m05_edges      L109_fibre_all_m05 + edges held on all 14 lofts
  L110_fibre_all_m02_edges      L109_fibre_all_m02 + the same
  L110_fibre_parcus_m05_edges   L109_fibre_parcus_m05 + the P-arcus lofts' edges held
  L110_fibre_usl_m05_edges      L109_fibre_usl_m05 + the USL lofts' edges held
usage: py -3.10 build_batch110.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch109 import Model109, fibre, ALL, BASE, RUNS, emit  # noqa: E402

WANT = set(sys.argv[1:])
BUILDS = [('L110_fibre_all_m05_edges', BASE, fibre(ALL, 0.05, edges=True)),
          ('L110_fibre_all_m02_edges', BASE, fibre(ALL, 0.02, edges=True)),
          ('L110_fibre_parcus_m05_edges', BASE, fibre(('parcus',), 0.05, edges=True)),
          ('L110_fibre_usl_m05_edges', BASE, fibre(('usl',), 0.05, edges=True))]

if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        if os.path.exists(os.path.join(RUNS, name)):
            print(f'{name} exists; not overwriting')
            continue
        mdl = Model109(os.path.join(RUNS, base, base + '.feb'))
        fn(mdl)
        emit(name, mdl)
