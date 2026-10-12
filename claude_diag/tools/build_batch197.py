"""Batch 197 (2026-10-07 night; tube line, secondary to the band fit). L192_botfull_distpress_cdamp1000_P1 (L187: the user's
equal PM-arch springs; + curve pressure on the distal half of the curves; curve damping C 1000) ran to t 1 with the smallest
fold of any tube run (3-10 of 828 nodes inside) at -5.5 / -9.7 / -37.3: Ba above the hymen, Bp 14 mm short, Bp - Ba -4.2.
The posterior support is the Bp lever on the band line (x1 -> x0.5: Bp +7 to +9, Ba +1), so one change on it, P1:
  L197_distpress_parcus04_P1   Parcus_conn x 0.4 (applied to the P1 file, i.e. healthy posterior x0.4 under P1's x0.15)
usage: py -3.10 build_batch197.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch145 import healthy  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L192_botfull_distpress_cdamp1000_P1'
NAME = 'L197_distpress_parcus04_P1'

if __name__ == '__main__':
    assert not os.path.exists(os.path.join(RUNS, NAME))
    m = Model7(os.path.join(RUNS, SRC, SRC + '.feb'))
    m.log.append(f'{NAME}: {SRC} (already P1) + healthy posterior (Parcus_conn) x0.4 (tools/build_batch197.py)')
    healthy(m, {'parcus': 0.4})
    emit(NAME, m)
