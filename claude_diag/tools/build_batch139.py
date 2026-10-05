"""Batch 139 (2026-10-01, the user's choice "2": rerun L138_tube_thin275_narrow3_flare_pmin_wrapP (the pressure on the curved
edges; stopped at t 0.25 crawling, steps 1e-4..1e-3, 344 'max iterations / reformations reached') with solver-only changes,
one each, checked against L138's states up to t 0.25 so the answer does not move:
  L139_wrapP_optiter25   time_stepper opt_iter 10 -> 25 (the auto-stepper stops cutting steps that take 11-25 iterations)
  L139_wrapP_dtol01      solver dtol 0.001 -> 0.01 (looser displacement norm; earlier lines: < 1 mm change)
usage: py -3.10 build_batch139.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L138_tube_thin275_narrow3_flare_pmin_wrapP'
V = {'L139_wrapP_optiter25': ('opt_iter', '10', '25'), 'L139_wrapP_dtol01': ('dtol', '0.001', '0.01')}
t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
for name, (tag, old, new) in V.items():
    d = os.path.join(RUNS, name)
    assert not os.path.exists(d), name
    a, b = f'<{tag}>{old}</{tag}>', f'<{tag}>{new}</{tag}>'
    assert t.count(a) == 1, (tag, t.count(a))
    os.makedirs(d)
    open(os.path.join(d, name + '.feb'), 'w', encoding='utf-8', newline='').write(t.replace(a, b))
    open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8').write(f'{name}: {SRC} + solver only: {tag} {old} -> {new}\n')
    print(name, tag, old, '->', new)
