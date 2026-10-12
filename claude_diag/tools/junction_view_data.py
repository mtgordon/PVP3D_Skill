"""Export the tube's distal junction (curves' PVW-side ends / PVW / PeB, near the hymen) for the 3D viewer (2026-10-10, the user:
"show me the junction picture first"). Canal frame as apex_view_data.py (x lateral, y AP, z depth from the PM plane). Writes the
outer faces of every part with a node in the box (depth < ZMAX), the curve hexes that sit INSIDE a PVW / PeB hex (centroid inside:
'curve_overlap', the overlap the point ties live in), and the canal_wrap_embed tied nodes as points. Optionally --at RUN T: the
deformed shape at the last converged state at or before T (from the run's .xplt).
usage: py -3.10 junction_view_data.py RUN OUT.json [--at RUN T]"""
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apex_close import HF  # noqa: E402
from build_batch202 import plane  # noqa: E402
from build_batch143 import natural  # noqa: E402
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
DOMS = {'_PickedSet347': 'AVW', '_PickedSet64': 'PVW', '_PickedSet66': 'PeB',
        'canal_wrap_hex_avw': 'curve_avw', 'canal_wrap_hex_pvw': 'curve_pvw'}
ZMAX = 14.0


def frame(X, E):
    c0, nrm = plane()
    cx = [v for n, c in E.values() if n == '_PickedSet346' for v in c]
    nn = nrm if np.dot(np.mean([X[v] for v in cx], 0) - c0, nrm) < 0 else -nrm
    e3 = -nn
    e1 = np.array([1.0, 0, 0]); e1 -= (e1 @ e3) * e3; e1 /= np.linalg.norm(e1)
    return c0, np.array([e1, np.cross(e3, e1), e3]).T


def load(run):
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {int(el.get('id')): (e.get('name'), [int(v) for v in el.text.split(',')])
         for e in mesh.findall('Elements') for el in e}
    return t, X, E


def deformed(X, run, tw):
    sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
    from xplt_reader import Xplt, converged_state_indices
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = [i for i in converged_state_indices(x, os.path.join(RUNS, run, run + '.log')) if x.states[i][0] <= tw + 1e-9][-1]
    U = x.var(idx, 'displacement')
    ids = getattr(x, 'node_ids', None)
    pos = {int(n): i for i, n in enumerate(ids)} if ids is not None else None
    return {k: p + U[pos[k] if pos else k - 1] for k, p in X.items()}, x.states[idx][0]


def overlap(X, E):
    """curve_pvw / curve_avw hexes whose centroid lies inside a PVW / PeB hex"""
    hosts = [(k, np.array([X[v] for v in c])) for k, (n, c) in E.items() if n in ('_PickedSet64', '_PickedSet66')]
    HC = np.array([P.mean(0) for _, P in hosts])
    out = set()
    for k, (n, c) in E.items():
        if not n.startswith('canal_wrap_hex'):
            continue
        p = np.mean([X[v] for v in c], 0)
        for j in np.argsort(np.linalg.norm(HC - p, axis=1))[:8]:
            if np.abs(natural(hosts[j][1], p)).max() <= 1.0:
                out.add(k)
                break
    return out


def main(run, out, at=None):
    t, X, E = load(run)
    c0, B = frame(X, E)
    over = overlap(X, E)
    tag = run
    if at:
        X, tt = deformed(X, at[0], float(at[1]))
        tag = f'{at[0]} at t {tt:.3f}'
    L = {k: (p - c0) @ B for k, p in X.items()}
    blk = re.search(r'<constraint name="(?:canal_wrap_embed|junction_ties)".*?</constraint>', t, re.S)
    tied = sorted({int(m) for m in re.findall(r'<linear_constraint>\s*<node id="(\d+)" bc="x">', blk.group(0))}) if blk else []
    Lr = {k: (p - c0) @ B for k, p in load(run)[1].items()}       # the box is chosen on the undeformed geometry
    parts, used = {}, set()
    for dom, lab in DOMS.items():
        els = [(k, c) for k, (n, c) in E.items() if n == dom]
        cnt = Counter(frozenset(c[i] for i in fc) for _, c in els for fc in HF)
        for k, c in els:
            if min(Lr[v][2] for v in c) > ZMAX:
                continue
            key = 'curve_overlap' if k in over else lab
            for fc in HF:
                f = [c[i] for i in fc]
                if cnt[frozenset(f)] == 1 or key == 'curve_overlap':
                    parts.setdefault(key, []).append(f)
                    used.update(f)
    used.update(tied)
    ids = sorted(used)
    idx = {v: i for i, v in enumerate(ids)}
    data = {'run': tag, 'nodes': [[round(float(c), 3) for c in L[v]] for v in ids], 'node_ids': ids,
            'parts': {k: [[idx[v] for v in f] for f in q] for k, q in parts.items()},
            'points': {'tied': [idx[v] for v in tied]}}
    json.dump(data, open(out, 'w'))
    print(out, tag, len(ids), 'nodes', {k: len(q) for k, q in parts.items()}, 'tied', len(tied), 'overlap hexes', len(over))


if __name__ == '__main__':
    a = sys.argv[1:]
    at = a[a.index('--at') + 1:a.index('--at') + 3] if '--at' in a else None
    main(a[0], a[1], at)
