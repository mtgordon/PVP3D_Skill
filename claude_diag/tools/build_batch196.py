"""Batch 196 (2026-10-07 night, the unattended fit). On the no-PVW-band base (L183: 6 mm band on the AVW + cervix) with
posterior x0.35, apical x6 gave P1 +1.6 / +3.8 and x8 -1.9 / +0.4: Bp - Ba stays ~+2.2; the paper needs >= ~+4. The band's
width on the AVW is the lever that holds Ba without the PVW (3 -> 6 mm on all walls: Ba -5.1, Bp -2.3), so widen it on the
AVW + cervix only, the PVW still band-free (seam-band changes only; the band now covers most of the AVW, said plainly):
  L196_band9av_0pv_pmd_sph025   = build_batch175.make with band(9, 3) on the AVW + cervix, none on the PVW
usage: py -3.10 build_batch196.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch175 import make  # noqa: E402

if __name__ == '__main__':
    make('L196_band9av_0pv_pmd_sph025', [(9, 3, ('Vagina_AVW', 'Vagina_Cervix'))], 'tools/build_batch196.py')
