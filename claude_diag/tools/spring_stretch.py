"""Elongation of each discrete-spring family at a plot state against the end of its force table (2026-10-07: are the
apical CL / USL springs past the last point of their source table, where the force stays constant?).
usage: py -3.10 spring_stretch.py RUN [--t 1] [--fam CL,USL]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt  # noqa: E402
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
args = sys.argv[1:]
tw, fams = 1.0, None
if '--t' in args:
    i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]
if '--fam' in args:
    i = args.index('--fam'); fams = args[i + 1].split(','); args = args[:i] + args[i + 2:]
run = args[0]
t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
X0 = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(',')])
      for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t)}
sets = {m.group(1): [tuple(int(v) for v in e.split(',')) for e in re.findall(r'<delem>(\d+,\d+)</delem>', m.group(2))]
        for m in re.finditer(r'<DiscreteSet name="([^"]+)">(.*?)</DiscreteSet>', t, re.S)}
# discrete materials in list order; <discrete dmat="k" discrete_set="name"> maps sets to them (by list position)
dm = re.findall(r'<discrete_material id="\d+" name="([^"]+)" type="([^"]+)">(.*?)</discrete_material>', t, re.S)
links = re.findall(r'<discrete dmat="(\d+)" discrete_set="([^"]+)"', t)
x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
si = min(range(len(x.states)), key=lambda i: abs(x.states[i][0] - tw))
u = np.asarray(x.var(si, 'displacement'), dtype=float)
idx = {int(n): i for i, n in enumerate(x.node_ids)}
print(f'{run} at t {x.states[si][0]:.3f}')
print('| set | material | n | elongation median / max (mm) | table end x (mm) | past the end | scale |')
print('|---|---|---|---|---|---|---|')
for k, sname in links:
    if fams and not any(sname.startswith(f) for f in fams):
        continue
    mname, mtype, body = dm[int(k) - 1]
    pts = [tuple(float(v) for v in p.split(',')) for p in re.findall(r'<pt>([^<]+)</pt>', body)]
    sc = re.search(r'<scale>([^<]+)</scale>', body)
    xend = max(p[0] for p in pts) if pts else float('nan')
    el = []
    for a, b in sets[sname]:
        L0 = np.linalg.norm(X0[b] - X0[a])
        L = np.linalg.norm(X0[b] + u[idx[b]] - X0[a] - u[idx[a]])
        el.append(L - L0)
    el = np.array(el)
    print(f'| {sname} | {mname} | {len(el)} | {np.median(el):.1f} / {el.max():.1f} | {xend:g} | {(el > xend).sum()} | '
          f'{sc.group(1) if sc else "-"} |')
