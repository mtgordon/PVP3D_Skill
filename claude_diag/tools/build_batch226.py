"""Batch 226 (2026-10-08 night, unattended; the user: "see if you can speed up the runs some without causing other issues").
Reference: L224_sl01_cl30parcus062_src_P1 (CL x30, Parcus x0.62, lips x0.1; +0.5 / +4.9 / -38.7; 2 h 26 min, 215 failed attempts,
most with negative Jacobians, 184 of them after t 0.8). One solver / output change per variant (the model unchanged):
  E1 diag      time_steps 100 -> 90 (to t 0.9), run through the console with output_negative_jacobians (which elements invert)
  E2 dtmax01   time_stepper dtmax 0.005 -> 0.01
  E3 stride2   plot_stride 1 -> 2 (half the plot output)
usage: py -3.10 build_batch226.py MODE"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L224_sl01_cl30parcus062_src_P1'
EDITS = {'diag': (r'<time_steps>100</time_steps>', '<time_steps>90</time_steps>',
                  'diagnostic: time_steps 100 -> 90 (to t 0.9), run with output_negative_jacobians'),
         'dtmax01': (r'<dtmax>0\.005</dtmax>', '<dtmax>0.01</dtmax>', 'solver only: dtmax 0.005 -> 0.01'),
         'stride2': (r'<plot_stride>1</plot_stride>', '<plot_stride>2</plot_stride>', 'output only: plot_stride 1 -> 2')}

if __name__ == '__main__':
    mode = sys.argv[1]
    pat, rep, msg = EDITS[mode]
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    t, n = re.subn(pat, rep, t, count=1)
    assert n == 1, mode
    write(f'L226_cl30p062_{mode}_P1', t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch226.py)\n{msg}\n')
