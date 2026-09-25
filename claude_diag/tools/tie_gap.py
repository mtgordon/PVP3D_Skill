"""Tie gaps of the tied-node-on-facet contacts at t = 1 (or the last converged state).

usage: py -3.10 tie_gap.py RUN [RUN ...]
Each tnof contact ties its primary nodes (tissue, _PickedSet surfaces) onto the secondary facets (the loft). As meshed,
every primary node starts on a loft node, so the gap is the distance between that pair at the reported state
(an upper bound on the distance to the loft surface). Also the tissue displacement there, for scale.
"""
import os
import sys

import numpy as np

import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR
from xplt import Xplt
from febmodel import Feb


def surface_nodes(feb, name):
    s = next(s for s in feb.root.find('Mesh').findall('Surface') if s.get('name') == name)
    return sorted({int(v) for f in s for v in f.text.split(',')})


def main(runs):
    for r in runs:
        f = Feb(os.path.join(RUNS_DIR, r, r + '.feb'))
        x = Xplt(os.path.join(RUNS_DIR, r, r + '.xplt'))
        idx = {int(nd): i for i, nd in enumerate(x.node_ids)}
        conv = [i for i, s in enumerate(x.states) if s[1] == 0]
        times = np.array([x.states[i][0] for i in conv])
        k = int(np.argmin(abs(times - 1.0))) if times[-1] >= 0.99 else len(times) - 1
        u = x.var(conv[k], 'displacement')
        pairs = {p.get('name'): (p.find('primary').text, p.find('secondary').text)
                 for p in f.root.find('Mesh').findall('SurfacePair')}
        for c in f.root.find('Contact'):
            if c.get('type') != 'tied-node-on-facet':
                continue
            prim, sec = pairs[c.get('surface_pair')]
            pn, sn = surface_nodes(f, prim), surface_nodes(f, sec)
            S = np.array([f.nodes[n] for n in sn])
            gaps, disp = [], []
            for n in pn:
                d = np.linalg.norm(S - f.nodes[n], axis=1)
                j = int(np.argmin(d))
                if d[j] > 1e-6:
                    continue
                gaps.append(np.linalg.norm(u[idx[n]] - u[idx[sn[j]]]))
                disp.append(np.linalg.norm(u[idx[n]]))
            g, dd = np.array(gaps), np.array(disp)
            print(f'{r:10s} t {times[k]:.3f} {c.get("name"):38s} penalty {c.find("penalty").text:>6s}: {len(g)} pairs, '
                  f'gap median {np.median(g):.3f} / p90 {np.percentile(g, 90):.3f} / max {g.max():.3f} mm '
                  f'(tissue moved median {np.median(dd):.1f} mm)')


if __name__ == '__main__':
    main(sys.argv[1:])
