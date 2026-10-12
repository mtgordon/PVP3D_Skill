"""Batch 199 (2026-10-07 night, the unattended fit). The 9 mm band on the AVW + cervix (L196, 72 % of the AVW at c1 x3; PVW
band-free) + apical x6 + posterior x0.35 MEETS P1 (-0.8 / +4.1; Bp - Ba +4.9), where the 6 mm band gave +1.6 / +3.8. To see
how little band does it (the user prefers not to stiffen whole walls), two single changes on that base, P1:
  L199_band9x2av_0pv_pmd_sph025   the 9 mm AVW + cervix band at c1 x2 instead of x3 (softer, same coverage)
  L199_band75av_0pv_pmd_sph025    a 7.5 mm AVW + cervix band at c1 x3 (between the 6 and 9 mm)
usage: py -3.10 build_batch199.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch175 import make  # noqa: E402

AC = ('Vagina_AVW', 'Vagina_Cervix')
if __name__ == '__main__':
    make('L199_band9x2av_0pv_pmd_sph025', [(9, 2, AC)], 'tools/build_batch199.py')
    make('L199_band75av_0pv_pmd_sph025', [(7.5, 3, AC)], 'tools/build_batch199.py')
