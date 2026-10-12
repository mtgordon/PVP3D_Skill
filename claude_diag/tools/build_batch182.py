"""Batch 182 (2026-10-06 night, the unattended fit): the seam band per wall, to let Bp descend while Ba is held.
On the 6 mm band base (L175_band6_pmd_sph025) with apical x5 + posterior x0.4, the best P1 was +3.9 / -0.8 (Bp 4.8 mm short)
and the posterior support cannot go lower (x0.3 broke on both bands). The band's width is the Ba lever (3 -> 6 mm: Ba -5.1,
Bp -2.3), so keep the AVW + cervix band at 6 mm and give the PVW only the 3 mm band (seam-only changes, NOT IN SOURCE):
  L182_band6av_3pv_pmd_sph025   = build_batch175.make with band(6, 3) on the AVW + cervix and band(3, 3) on the PVW
usage: py -3.10 build_batch182.py   (then build_batch145.py ... --base L182_band6av_3pv_pmd_sph025)"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch175 import make  # noqa: E402

if __name__ == '__main__':
    make('L182_band6av_3pv_pmd_sph025', [(6, 3, ('Vagina_AVW', 'Vagina_Cervix')), (3, 3, ('Vagina_PVW',))],
         'tools/build_batch182.py')
