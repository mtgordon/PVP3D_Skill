"""Where is a failed FEBio run deforming, and where does it blow up?

Reads the run's .feb (for domain names) and .xplt (for results) and prints:
  1. max |u| per domain over the converged history (which part is running away?)
  2. relative volume J extremes per domain at the last converged state, with the
     element IDs (which element is closest to inverting?)
  3. if the plot file holds the failed iterations after the last converged state
     (FEBioStudio debug runs, or <plot_level>PLOT_MINOR_ITRS</plot_level>), the
     nodes whose increment blows up, grouped by domain -- this localizes a failure
     that the log only reports as "N negative jacobians detected".

usage: py -3 feb_postmortem.py model.feb model.xplt [model.log] [--iters 40]

The .log is needed only for plain `febio4 -i` runs that plotted minor iterations
(their states are not flagged converged; see xplt_reader.converged_state_indices).
Needs Python 3 + numpy; imports xplt_reader.py from this directory.
"""
import argparse
import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xplt_reader import Xplt, converged_state_indices  # noqa: E402


def domain_names(feb_path):
    root = ET.parse(feb_path).getroot()
    names = [d.get('name') for d in root.find('MeshDomains')]
    disc = root.find('Discrete')
    if disc is not None:
        names += ['discrete:' + d.get('discrete_set') for d in disc.findall('discrete')]
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('feb')
    ap.add_argument('xplt')
    ap.add_argument('log', nargs='?')
    ap.add_argument('--iters', type=int, default=40, help='failed iterations to inspect')
    ap.add_argument('--samples', type=int, default=8, help='converged states in the history table')
    a = ap.parse_args()

    x = Xplt(a.xplt)
    names = domain_names(a.feb)
    for i, d in enumerate(x.domains):
        d['name'] = names[i] if i < len(names) else f'region{i + 1}'
    solid_like = [i for i, d in enumerate(x.domains) if not d['name'].startswith('discrete:')]

    conv = converged_state_indices(x, a.log)
    if not conv:
        sys.exit('no converged states found (pass the .log for PLOT_MINOR_ITRS runs)')
    last = conv[-1]
    print(f'{len(x.states)} states, {len(conv)} converged; last converged t={x.states[last][0]:.6g}')

    # 1) displacement history per domain
    step = max(1, len(conv) // a.samples)
    sel = sorted(set(conv[::step] + [last]))
    hist = {}
    for s in sel:
        mag = np.linalg.norm(x.var(s, 'displacement'), axis=1)
        hist[s] = {i: float(mag[np.unique(x.domains[i]['conn'].ravel())].max()) for i in solid_like}
    print('\nmax |u| per domain at converged states (column = time):')
    print('%-30s' % 'domain' + ''.join('%10.4g' % x.states[s][0] for s in sel))
    for i in solid_like:
        print('%-30s' % x.domains[i]['name'][:30] + ''.join('%10.3f' % hist[s][i] for s in sel))

    # 2) J extremes at the last converged state
    J = x.var(last, 'relative volume')
    if J is not None:
        print('\nrelative volume J at last converged state (sorted by min J):')
        rows = []
        for i in solid_like:
            arr = J.get(i + 1)
            if arr is None or not len(arr):
                continue
            k, K = int(np.argmin(arr)), int(np.argmax(arr))
            rows.append((float(arr[k]), x.domains[i]['name'], int(x.domains[i]['eids'][k]),
                         float(arr[K]), int(x.domains[i]['eids'][K])))
        for r in sorted(rows)[:15]:
            print('  %-30s minJ %8.4f (elem %d)   maxJ %8.4f (elem %d)' % (r[1][:30], r[0], r[2], r[3], r[4]))

    # 3) failed iterations after the last converged state
    after = list(range(last + 1, min(len(x.states), last + 1 + a.iters)))
    if not after:
        print('\nno iteration states after the last converged one (plot with PLOT_MINOR_ITRS to see them)')
        return
    node_dom = defaultdict(set)
    for i in solid_like:
        for n in np.unique(x.domains[i]['conn'].ravel()):
            node_dom[int(n)].add(x.domains[i]['name'])
    u0 = x.var(last, 'displacement')
    print('\nfailed iterations after last converged state: largest increment and where the top-50 nodes live')
    for s in after:
        du = np.linalg.norm(x.var(s, 'displacement') - u0, axis=1)
        order = np.argsort(du)[::-1]
        cnt = defaultdict(int)
        for n in order[:50]:
            for dn in node_dom.get(int(n), {'(no element)'}):
                cnt[dn] += 1
        where = ', '.join(f'{k}:{v}' for k, v in sorted(cnt.items(), key=lambda kv: -kv[1])[:4])
        print('  state %5d t=%.6g  max|du|=%9.3g at node %d  [%s]' % (
            s, x.states[s][0], du[order[0]], x.node_ids[order[0]], where))


if __name__ == '__main__':
    main()
