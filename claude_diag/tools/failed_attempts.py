"""Localize every failed time-step attempt in a PLOT_MINOR_ITRS run (not just the ones after the
last converged state): for each failed attempt, the iteration with the largest increment relative
to the previous converged state, and which domains its top-50 nodes belong to.
usage: py -3.10 failed_attempts.py model.feb model.xplt model.log [--max 20]"""
import argparse, os, re, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt  # noqa: E402
import xml.etree.ElementTree as ET

ap = argparse.ArgumentParser(); ap.add_argument('feb'); ap.add_argument('xplt'); ap.add_argument('log')
ap.add_argument('--max', type=int, default=20); a = ap.parse_args()
x = Xplt(a.xplt)
names = [d.get('name') for d in ET.parse(a.feb).getroot().find('MeshDomains')]
node_dom = defaultdict(set)
for i, d in enumerate(x.domains):
    nm = names[i] if i < len(names) else f'region{i + 1}'
    for n in np.unique(d['conn'].ravel()):
        node_dom[int(n)].add(nm)
log = open(a.log, encoding='latin-1').read()
events = [(m.group(1), float(m.group(2))) for m in re.finditer(r'^------- (converged|failed to converge) at time : ([0-9.eE+-]+)', log, re.M)]
times = np.array([s[0] for s in x.states])
# walk states in order; a state whose time matches a converged time (and is the last with that time) is converged
conv_t = [t for k, t in events if k == 'converged']
conv_idx = []
for t in conv_t:
    js = np.where(np.abs(times - t) < 1e-6 * max(t, 1e-6))[0]
    if len(js):
        conv_idx.append(int(js.max()))
conv_idx = sorted(set(conv_idx))
failed_t = [t for k, t in events if k == 'failed to converge']
print(f'{len(x.states)} states, {len(conv_idx)} converged matched, {len(failed_t)} failed attempts')
for t in failed_t[:a.max]:
    js = np.where(np.abs(times - t) < 1e-6 * max(t, 1e-6))[0]
    if not len(js):
        print(f'  t={t:.6g}: no plotted iterations'); continue
    base = max([c for c in conv_idx if c < js.min()], default=None)
    if base is None:
        continue
    u0 = x.var(base, 'displacement')
    best = None
    for j in js:
        du = np.linalg.norm(x.var(int(j), 'displacement') - u0, axis=1)
        if best is None or du.max() > best[0]:
            best = (float(du.max()), int(j), du)
    mx, j, du = best
    order = np.argsort(du)[::-1][:50]
    cnt = defaultdict(int)
    for n in order:
        for dn in node_dom.get(int(n), {'(no element)'}):
            cnt[dn] += 1
    where = ', '.join(f'{k}:{v}' for k, v in sorted(cnt.items(), key=lambda kv: -kv[1])[:4])
    print(f'  failed t={t:.6g}: {len(js)} iters, max |du| {mx:9.3g} at node {x.node_ids[order[0]]} [{where}]')
