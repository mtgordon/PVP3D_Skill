"""Batch 37 (2026-09-25, user): the loft model with the LA pressure at ~1/3 and the other pressures at ~2x the source.

The user's values (rounded): Load-LA 0.005 MPa (5 kPa; the source 0.014, so ~0.36x; the lofts line had 0.00467 = 1/3)
and Load-AVW, Load-PVW, Load-PeB-top, Load-top 0.03 MPa (30 kPa, ~2.1x the source 0.014). NOT IN SOURCE (both).
  L27_lofts_LA5kPa_rest30kPa  L26_lofts_la3_pm (the lofts line, current version, with PM_Plane) + those pressures
usage: py -3.10 build_batch37.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L26_lofts_la3_pm', 'L26_lofts_la3_pm.feb')
NEW = {'Load-LA': 0.005, 'Load-AVW': 0.03, 'Load-PVW': 0.03, 'Load-PeB-top': 0.03, 'Load-top': 0.03}
WANT = set(sys.argv[1:])

name = 'L27_lofts_LA5kPa_rest30kPa'
if not WANT or name in WANT:
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(BASE)
    done = set()
    for sl in m.root.find('Loads'):
        if sl.get('name') in NEW:
            p = sl.find('pressure')
            old = p.text
            p.text = '%g' % NEW[sl.get('name')]
            done.add(sl.get('name'))
            m.log.append(f'NOT IN SOURCE (user): {sl.get("name")} pressure {old} -> {p.text} MPa '
                         f'({NEW[sl.get("name")] / 0.014:.3g}x the source 0.014)')
    assert done == set(NEW), done
    emit(name, m)
