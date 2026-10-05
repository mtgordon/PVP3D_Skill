"""Batch 147 (2026-10-02, the user: "if you know of things that lead to essentially the same result but run faster (maybe
removing some contacts), you can use those as intermediate cases to learn").
The fast copy of the paper-fit base: L144_seamspr_pm_allcontact without the NOT-IN-SOURCE wall contacts walls_LA and
walls_PM (their contacts, SurfacePairs and surfaces walls_LA_primary, walls_PM_primary, PM_contact removed). On the tube
line these cost 3-4x the run time in P1 / P2 and moved the answer 0.2 mm (healthy) to ~1-2 mm (P1 / P2). Everything else
(seam springs, seg_up 2, the PM structure, the 27 PM_conn-family springs on the inner arc, PVW_LA, SlidingElastic1) stays.
Used only to learn the sensitivities fast; the chosen set is confirmed on L144.
  L147_seamspr_pm_fast   built, test-read, then run
usage: py -3.10 build_batch147.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE, NAME = 'L144_seamspr_pm_allcontact', 'L147_seamspr_pm_fast'

if __name__ == '__main__':
    d = os.path.join(RUNS, NAME)
    assert not os.path.exists(d), f'{NAME} exists; not overwriting'
    t = open(os.path.join(RUNS, BASE, BASE + '.feb'), encoding='ISO-8859-1').read()
    n0 = len(t)
    for nm in ('walls_LA', 'walls_PM'):
        t, k = re.subn(r'[ \t]*<contact [^>]*name="%s"[^>]*>.*?</contact>\s*\n' % nm, '', t, flags=re.S)
        assert k == 1, (nm, k)
        t, k = re.subn(r'[ \t]*<SurfacePair name="%s">.*?</SurfacePair>\s*\n' % nm, '', t, flags=re.S)
        assert k == 1, (nm, k)
    for s in ('walls_LA_primary', 'walls_PM_primary', 'PM_contact'):
        t, k = re.subn(r'[ \t]*<Surface name="%s">.*?</Surface>\s*\n' % s, '', t, flags=re.S)
        assert k == 1, (s, k)
    for s in ('walls_LA', 'walls_PM', 'PM_contact'):
        assert s not in t, s
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='ISO-8859-1', newline='').write(t)
    msg = (f'base: {os.path.join(RUNS, BASE, BASE + ".feb")} (tools/build_batch147.py)\n{NAME}: {BASE} without the contacts '
           f'walls_LA and walls_PM (contacts, SurfacePairs, surfaces walls_LA_primary / walls_PM_primary / PM_contact '
           f'removed; {n0 - len(t)} characters); the fast copy for learning the paper-fit sensitivities\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
