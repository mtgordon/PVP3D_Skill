"""Solver effort per run from the log (Test A speed comparisons): converged steps, failed attempts, equilibrium
iterations, right-hand-side evaluations, stiffness reformations, equations, wall time, time in the linear solver, and
seconds per iteration. usage: py -3.10 speed_table.py RUN [RUN ...]"""
import os
import re
import sys

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def hms(s):
    h, m, sec = (int(v) for v in s.split(':'))
    return 3600 * h + 60 * m + sec


print('| run | t | steps | failed | iterations | RHS evals | reformations | equations | wall | linear solver | s / iteration |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for run in sys.argv[1:]:
    txt = open(os.path.join(RUNS, run, run + '.log'), encoding='latin-1').read()
    g = lambda pat: (re.findall(pat, txt) or [''])[-1]
    t = g(r'------- converged at time : ([0-9.eE+-]+)')
    steps = g(r'Number of time steps completed \.+ : (\d+)')
    its = g(r'Total number of equilibrium iterations \.+ : (\d+)')
    rhs = g(r'Total number of right hand evaluations \.+ : (\d+)')
    refs = g(r'Total number of stiffness reformations \.+ : (\d+)')
    neq = g(r'Nr of equations \.+ : (\d+)')
    wall = g(r'Total elapsed time \.+ : (\d+:\d+:\d+)')
    lin = g(r'time in linear solver \.+ : (\d+:\d+:\d+)')
    fails = len(re.findall(r'------- failed to converge', txt))
    spi = f'{hms(wall) / int(its):.3f}' if wall and its else ''
    print(f'| `{run}` | {t} | {steps} | {fails} | {its} | {rhs} | {refs} | {neq} | {wall} | {lin} | {spi} |')
