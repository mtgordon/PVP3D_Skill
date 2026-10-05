"""Batch 160 (2026-10-02, the user: "yes, try the seam springs x0.5 on C3").
C3 on the fast copy (L157_C3_parcus042_<case>: L149_seamspr_pm_outer_fast + healthy Parcus x0.42, then the case) with the canal
seam springs (constraint canal_seam_springs: 116 AVW / cervix edge nodes tied to the PVW edge, zero-length, penalty 0.01 N/mm
per dof, maxaug 0; NOT IN SOURCE) at half stiffness: penalty 0.01 -> 0.005. One change from L157_C3_parcus042_<case>.
  L160_C3_seamx05_healthy / _P1 / _P2
usage: py -3.10 build_batch160.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCALE = 0.5

if __name__ == '__main__':
    for case in ('healthy', 'P1', 'P2'):
        src, name = f'L157_C3_parcus042_{case}', f'L160_C3_seamx05_{case}'
        d = os.path.join(RUNS, name)
        assert not os.path.exists(d), f'{name} exists; not overwriting'
        p = os.path.join(RUNS, src, src + '.feb')
        t = open(p, encoding='ISO-8859-1').read()
        i = t.index('<constraint name="canal_seam_springs"')
        j = t.index('</constraint>', i)
        blk = t[i:j]
        old = re.search(r'<penalty>([^<]+)</penalty>', blk).group(1)
        new = '%g' % (float(old) * SCALE)
        blk2, k = re.subn(r'<penalty>[^<]+</penalty>', f'<penalty>{new}</penalty>', blk, count=1)
        assert k == 1
        t = t[:i] + blk2 + t[j:]
        os.makedirs(d)
        open(os.path.join(d, name + '.feb'), 'w', encoding='ISO-8859-1', newline='').write(t)
        prev = open(p + '.changes.txt', encoding='utf-8', errors='replace').read()
        with open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
            fo.write(f'{name} (tools/build_batch160.py) from {p}:\n{prev}\nNOT IN SOURCE (the seam springs themselves are not in '
                     f'the source): constraint canal_seam_springs penalty {old} -> {new} N/mm per dof (x{SCALE}; 116 nodes)\n')
        print('wrote', name, old, '->', new)
