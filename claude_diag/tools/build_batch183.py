"""Batch 183 (2026-10-06 night, the unattended fit): the PVW without a seam band.
L182 (6 mm band on the AVW + cervix, 3 mm on the PVW) moved P1 Bp +3.8 mm against the 6 mm band everywhere at t 0.79 with Ba
unchanged: the PVW band holds Bp up. One step further, the PVW faithful (no band; the seam is still contact-only):
  L183_band6av_0pv_pmd_sph025   = build_batch175.make with band(6, 3) on the AVW + cervix only (NOT IN SOURCE)
usage: py -3.10 build_batch183.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch175 import make  # noqa: E402

if __name__ == '__main__':
    make('L183_band6av_0pv_pmd_sph025', [(6, 3, ('Vagina_AVW', 'Vagina_Cervix'))], 'tools/build_batch183.py')
