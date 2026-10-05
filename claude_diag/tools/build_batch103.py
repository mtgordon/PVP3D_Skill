"""Batch 103 (2026-09-28 ~16:56, the user away): can a solver tolerance carry the c1 x1.5 walls through?

L99_newline_vwyeoh_rhoi0_c1x2 reached t = 1 (5 failed); L99_newline_vwyeoh_rhoi0_c1x15 passed the faithful walls' crawl
(t 0.78) but crawls at t ~0.80-0.82 (~0.005 per 5 min). On the faithful walls the displacement tolerance dtol 0.001 ->
0.01 went furthest of the settings tried (L86_newline_vwyeoh_rhoi0_dtol01: t 0.808 vs 0.781) with the same solution
(within 0.2 mm at t 0.76): the rattling seam nodes' increments keep the displacement norm from meeting 0.001.
  L103_newline_vwyeoh_rhoi0_c1x15_dtol01   L99_newline_vwyeoh_rhoi0_c1x15 + solver dtol 0.001 -> 0.01 (a tolerance only)
17:04: c1 x2 and x2.5 reached t = 1 on the springs line. Does c1 x2 hold on the lofts line (the faithful walls crawled
there from t ~0.52: L87_lofts_newline_vwyeoh_rhoi0 = the lofts line L87_lofts_newline_rhoi0 + the Yeoh walls)?
  L104_lofts_newline_vwyeoh_rhoi0_c1x2     L87_lofts_newline_vwyeoh_rhoi0 + the three walls' Yeoh c1 x2 (NOT FAITHFUL)
17:52: c1 x1.5 + dtol 0.01 ran fast past t 0.92 (x1.5 alone slowed at 0.97-0.98); the same tolerance on x1.25:
  L103_newline_vwyeoh_rhoi0_c1x125_dtol01  L101_newline_vwyeoh_rhoi0_c1x125 + solver dtol 0.001 -> 0.01 (a tolerance only)
18:05: c1 x1.5 + dtol 0.01 reached t = 1 on the springs line (dtol: the same solution). The working set carried to the
lofts line (as L87_lofts_newline_vwyeoh_rhoi0 carried the faithful walls; not a one-change variant):
  L105_lofts_newline_vwyeoh_rhoi0_c1x15_dtol01  L87_lofts_newline_vwyeoh_rhoi0 + the walls' c1 x1.5 + solver dtol 0.01
18:22: c1 x2 on the lofts line crawls at t ~0.60 (the distal AVW seam nodes and the PVW_LA uterosacral ends). The bracket:
  L104_lofts_newline_vwyeoh_rhoi0_c1x3     L87_lofts_newline_vwyeoh_rhoi0 + the three walls' Yeoh c1 x3 (NOT FAITHFUL)
18:44: the seam band (c1 x3 within 5 mm of the seam, the rest faithful) converged on the springs line; on the lofts line:
  L106_lofts_newline_vwyeoh_rhoi0_band5x3  L87_lofts_newline_vwyeoh_rhoi0 + build_batch100.band(5, 3) (NOT IN SOURCE)
usage: py -3.10 build_batch103.py [NAME ...]   (an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch99 import c1x
from build_batch100 import band

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def ctrl(note, **kw):
    def fn(m):
        m.set_control(**kw)
        m.log.append(note)
    return fn


BUILDS = (('L103_newline_vwyeoh_rhoi0_c1x15_dtol01', 'L99_newline_vwyeoh_rhoi0_c1x15',
           ctrl('solver only: dtol 0.001 -> 0.01 (convergence tolerance); the model is unchanged', dtol=0.01)),
          ('L104_lofts_newline_vwyeoh_rhoi0_c1x2', 'L87_lofts_newline_vwyeoh_rhoi0', c1x(2)),
          ('L103_newline_vwyeoh_rhoi0_c1x125_dtol01', 'L101_newline_vwyeoh_rhoi0_c1x125',
           ctrl('solver only: dtol 0.001 -> 0.01 (convergence tolerance); the model is unchanged', dtol=0.01)),
          ('L104_lofts_newline_vwyeoh_rhoi0_c1x3', 'L87_lofts_newline_vwyeoh_rhoi0', c1x(3)),
          ('L106_lofts_newline_vwyeoh_rhoi0_band5x3', 'L87_lofts_newline_vwyeoh_rhoi0', band(5, 3)),
          ('L105_lofts_newline_vwyeoh_rhoi0_c1x15_dtol01', 'L87_lofts_newline_vwyeoh_rhoi0',
           lambda m: (c1x(1.5)(m), ctrl('solver only: dtol 0.001 -> 0.01 (convergence tolerance)', dtol=0.01)(m))))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
