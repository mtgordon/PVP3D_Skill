"""Batch 137 (2026-10-01, the user: "I don't like how narrow the vaginal wall is now. Instead of bringing it in 6 on each
side, bring it in only 3mm on each side and let me see it"). L136_tube_thin275_flare_pmin (build_batch136) with the
wall narrowed 3 mm on each side instead of 6 (wall_reshape.NARROW = 3); everything else the same: 2.75 mm walls with
their outer surfaces kept, the PVW flaring back to full width toward the perineal body, the curved edges 2.75 mm thick in
AVW / PVW halves, PVW_LA + walls_LA contacts, seg_up 2, the PM as a structure, PM_conn on the PM's inner arc (mean-scaled).
  L137_tube_thin275_narrow3_flare_pmin   built only
usage: py -3.10 build_batch137.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wall_reshape  # noqa: E402
import build_batch136 as b136  # noqa: E402

if __name__ == '__main__':
    wall_reshape.NARROW = 3.0
    b136.run_main(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_batch136.py'),
                  {'NAME': 'L137_tube_thin275_narrow3_flare_pmin', 'PRE': 'L137tmp_pre', 'MID': 'L137tmp_mid'})
