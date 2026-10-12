"""Batch 211 (2026-10-07 night, unattended): the .inp AS WRITTEN. The source is not a healthy model: its levator is
LA_Yamada50% and its loads match a published figure of the Abaqus model made with 50 % levator impairment (SESSION_SUMMARY,
Attempt 9). With the lips the source supports now run (L209_src_sl01_healthy reached t 1), so the literal replica can be run:
L209_src_sl01_healthy (source connector values, healthy levator, Load-LA 1/3) with
  1. the levator at 50 % of the healthy Yamada law (LA_Yamada50%, as the .inp; build_batch107.la_healthy_pct);
  2. Load-LA at the source's 0.014 MPa.
Kept from the line (standing rules, NOT IN SOURCE): PM structure + PM_conn on the outer arc, sphincter side x0.25, wall /
PM damping, soft seam lips x0.1, arcus chain mass x10.
usage: py -3.10 build_batch211.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch107 import la_healthy_pct  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L209_src_sl01_healthy'
NAME = 'L211_inp_as_written'

if __name__ == '__main__':
    tmp = NAME + '_tmp'
    m = Model7(os.path.join(RUNS, SRC, SRC + '.feb'))
    la_healthy_pct(m, 50)
    emit(tmp, m)
    p = os.path.join(RUNS, tmp, tmp + '.feb')
    t = open(p, encoding='utf-8').read()
    notes = open(p + '.changes.txt', encoding='utf-8').read()
    t, n = re.subn(r'(<surface_load name="Load-LA" surface="Load-LA_surf" type="pressure">\s*<pressure lc="1">)'
                   r'[^<]+(</pressure>)', r'\g<1>0.014\2', t)
    assert n == 1
    write(NAME, t, notes + 'as the source: Load-LA 0.014 MPa (was 1/3)\n(tools/build_batch211.py: the .inp as written)\n')
    import shutil
    shutil.rmtree(os.path.join(RUNS, tmp))
