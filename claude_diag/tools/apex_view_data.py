"""Export the tube's apical corners for the 3D viewer (2026-10-09, the user: "Render first, I'll react").
Writes the outer (boundary) faces of the walls, the curves and any apex parts within a box around both lateral apical corners,
plus the lumen-surface edges (SlidingElastic1Secondary = cervix / AVW lumen, SlidingElastic1Primary = PVW lumen), in the
canal frame (x lateral, y AP, z depth from the PM plane, build_batch202.plane()), as JSON.
usage: py -3.10 apex_view_data.py RUN OUT.json"""
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
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
DOMS = {'_PickedSet346': 'cervix', '_PickedSet347': 'AVW', '_PickedSet64': 'PVW',
        'canal_wrap_hex_avw': 'curve_avw', 'canal_wrap_hex_pvw': 'curve_pvw',
        'canal_wrap_hex_avw_apex': 'apex_avw', 'canal_wrap_hex_pvw_apex': 'apex_pvw'}
ZMIN, XMIN = 72.0, 6.0          # the box: depth > ZMIN, |lateral| > XMIN


def main(run, out):
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    c0, nrm = plane()
    cxn = [int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet346' for el in e for v in el.text.split(',')]
    nn = nrm if np.dot(np.mean([X[v] for v in cxn], 0) - c0, nrm) < 0 else -nrm
    e3 = -nn
    e1 = np.array([1.0, 0, 0]); e1 -= (e1 @ e3) * e3; e1 /= np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    B = np.array([e1, e2, e3]).T
    L = {k: (p - c0) @ B for k, p in X.items()}
    used, parts = {}, {}
    for e in mesh.findall('Elements'):
        name = e.get('name')
        if name not in DOMS:
            continue
        els = [[int(v) for v in el.text.split(',')] for el in e]
        cnt = Counter(frozenset(el[i] for i in fc) for el in els for fc in HF)
        quads = []
        for el in els:
            for fc in HF:
                f = [el[i] for i in fc]
                if cnt[frozenset(f)] != 1:
                    continue
                P = np.array([L[v] for v in f])
                if P[:, 2].max() < ZMIN or np.abs(P[:, 0]).max() < XMIN:
                    continue
                quads.append(f)
                for v in f:
                    used[v] = True
        parts[DOMS[name]] = quads
    edges = {}
    for sname, lab in (('SlidingElastic1Secondary', 'edge_top'), ('SlidingElastic1Primary', 'edge_bottom')):
        b = re.search(r'<Surface name="%s">(.*?)</Surface>' % sname, t, re.S).group(1)
        F = [[int(v) for v in x.split(',')] for x in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', b)]
        ec = Counter(tuple(sorted((f[k], f[(k + 1) % len(f)]))) for f in F for k in range(len(f)))
        segs = [e for e, c in ec.items() if c == 1 and max(L[e[0]][2], L[e[1]][2]) > ZMIN]
        for e in segs:
            for v in e:
                used[v] = True
        edges[lab] = [list(e) for e in segs]
    ids = sorted(used)
    idx = {v: i for i, v in enumerate(ids)}
    data = {'run': run, 'nodes': [[round(float(c), 3) for c in L[v]] for v in ids], 'node_ids': ids,
            'parts': {k: [[idx[v] for v in f] for f in q] for k, q in parts.items()},
            'edges': {k: [[idx[v] for v in e] for e in s] for k, s in edges.items()}}
    json.dump(data, open(out, 'w'))
    print(out, len(ids), 'nodes', {k: len(q) for k, q in parts.items()}, {k: len(s) for k, s in edges.items()})


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
