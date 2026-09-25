"""Post-mortem of a failed FEBio run from its .xplt: where is the model deforming/blowing up?

usage: postmortem.py model.feb model.xplt
"""
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

from xplt import Xplt

feb_path, xplt_path = sys.argv[1], sys.argv[2]
root = ET.parse(feb_path).getroot()
dom_names = [d.get('name') for d in root.find('MeshDomains')]
disc_names = [d.get('discrete_set') for d in root.find('Discrete').findall('discrete')]

x = Xplt(xplt_path)
names = []
for i, d in enumerate(x.domains):
    if i < len(dom_names):
        names.append(dom_names[i])
    else:
        j = i - len(dom_names)
        names.append('discrete:' + (disc_names[j] if j < len(disc_names) else str(j)))
for d, n in zip(x.domains, names):
    d['name'] = n

# node -> domains
node_dom = defaultdict(set)
id2idx = {int(n): i for i, n in enumerate(x.node_ids)}
for d in x.domains:
    if d['name'].startswith('discrete:stab_'):
        continue
    for c in d['conn']:
        for n in c:
            node_dom[int(n)].add(d['name'])

conv = [i for i, s in enumerate(x.states) if s[1] == 0]
print(f'{len(x.states)} states, {len(conv)} converged; last converged t={x.states[conv[-1]][0]:.6f}')

# 1) displacement history per domain (max |u|) over converged states
def per_dom_max(u):
    mag = np.linalg.norm(u, axis=1)
    out = {}
    for d in x.domains:
        if d['name'].startswith('discrete:'):
            continue
        idx = np.unique(d['conn'].ravel())
        # plot-file connectivity is 0-based node index
        out[d['name']] = float(mag[idx].max())
    return out

sel = conv[::max(1, len(conv) // 8)] + [conv[-1]]
hist = {}
for i in sel:
    u = x.var(i, 'displacement')
    hist[i] = per_dom_max(u)
print('\nmax |u| per domain at converged states (t):')
cols = sel
print('%-28s' % 'domain' + ''.join('%9.4f' % x.states[i][0] for i in cols))
for n in hist[cols[0]]:
    print('%-28s' % n[:28] + ''.join('%9.3f' % hist[i][n] for i in cols))

# 2) relative volume extremes at the last converged state
i_last = conv[-1]
J = x.var(i_last, 'relative volume')
print('\nrelative volume (J) extremes at last converged state:')
for di, d in enumerate(x.domains):
    if d['name'].startswith('discrete:'):
        continue
    arr = J.get(di + 1) if J is not None else None
    if arr is None:
        continue
    k = int(np.argmin(arr)); K = int(np.argmax(arr))
    print('  %-28s minJ=%8.4f (elem %d)  maxJ=%8.4f (elem %d)' % (d['name'][:28], arr[k], d['eids'][k], arr[K], d['eids'][K]))

# 3) failed iterations after last converged: where does the increment blow up?
u0 = x.var(i_last, 'displacement')
after = list(range(i_last + 1, min(len(x.states), i_last + 40)))
print('\nfailed-iteration increments after last converged state:')
for i in after:
    u = x.var(i, 'displacement')
    du = np.linalg.norm(u - u0, axis=1)
    top = np.argsort(du)[::-1][:6]
    doms = Counter = defaultdict(int)
    for n in np.argsort(du)[::-1][:50]:
        for dn in node_dom.get(int(n), {'(no elem)'}):
            doms[dn] += 1
    dsum = ', '.join(f'{k}:{v}' for k, v in sorted(doms.items(), key=lambda kv: -kv[1])[:4])
    print('  state %d t=%.6f  max|du|=%9.3f at node %d  top50 in [%s]' % (
        i, x.states[i][0], du[top[0]], x.node_ids[top[0]], dsum))
