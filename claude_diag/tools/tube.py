"""The vaginal canal as a tube (2026-09-29, the user: "create a model that actually connects the AVW and PVW so they form a
tube instead of two separate sheets. Wrap it around similar to how the cervix wraps around to the pvw"; both a tight fold
and a rounded wrap, compared).
Geometry helpers: the canal's lateral seam (tools/build_batch90.seam_nodes: the AVW / cervix lumen-edge nodes along both
sides), each wall's rim (its boundary faces perpendicular to the wall, i.e. the side faces) and, for every lumen-edge node,
its column of rim nodes through the wall thickness (lumen side first).
usage: py -3.10 tube.py [RUN]   (prints the seam and column structure)"""
import os
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from febmodel import Feb  # noqa: E402
from paths import RUNS_DIR  # noqa: E402

HEXF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
HEXE = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
AVW, PVW, CX, PEB = '_PickedSet347', '_PickedSet64', '_PickedSet346', '_PickedSet66'


def boundary_faces(f, doms):
    cnt, owner = Counter(), {}
    for dom in doms:
        for eid, c in f.elem_blocks[dom][1].items():
            for fc in HEXF:
                k = tuple(sorted(c[i] for i in fc))
                cnt[k] += 1
                owner[k] = [c[i] for i in fc]
    return [owner[k] for k, v in cnt.items() if v == 1]


def adjacency(f, doms):
    adj = defaultdict(set)
    for dom in doms:
        for c in f.elem_blocks[dom][1].values():
            for i, j in HEXE:
                adj[c[i]].add(c[j])
                adj[c[j]].add(c[i])
    return adj


def lumen_normal_field(f, surf):
    """(node -> unit normal of the lumen surface, averaged over its facets) for a named surface."""
    acc = defaultdict(lambda: np.zeros(3))
    for t, c in f.surfaces[surf]:
        P = np.array([f.nodes[n] for n in c])
        nn = np.cross(P[2] - P[0], P[3] - P[1]) if len(c) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
        for n in c:
            acc[n] += nn
    return {n: v / np.linalg.norm(v) for n, v in acc.items()}


def columns(f, doms, surf, edge_nodes):
    """For each lumen-edge node: the chain of nodes through the wall (lumen side first), walking hex edges on the
    wall's rim (boundary nodes) away from the lumen surface: at each step the neighbour with the largest gain along the
    local lumen normal (pointing into the wall), until no neighbour gains more than 0.3 mm."""
    adj = adjacency(f, doms)
    bf = boundary_faces(f, doms)
    bnodes = {n for fc in bf for n in fc}
    nrm = lumen_normal_field(f, surf)
    out = {}
    for v in edge_nodes:
        n0 = nrm.get(v)
        if n0 is None:
            continue
        # the lumen normal points out of the lumen surface's facets; into the wall = the side with the wall's nodes
        wall_side = np.mean([f.nodes[w] - f.nodes[v] for w in adj[v]], axis=0)
        d = n0 if wall_side @ n0 > 0 else -n0
        col, cur = [v], v
        while True:
            cand = [(float((f.nodes[w] - f.nodes[cur]) @ d), w) for w in adj[cur] if w in bnodes and w not in col]
            cand = [c for c in cand if c[0] > 0.3]
            if not cand:
                break
            gain, nxt = max(cand)
            col.append(nxt)
            cur = nxt
        out[v] = col
    return out


def main(run='L91_newline_vwyeoh_rhoi0_seamspr001'):
    from variants7 import Model7
    from build_batch90 import seam_nodes, _bedges
    path = os.path.join(RUNS_DIR, run, run + '.feb')
    f = Feb(path)
    seam, fP, X = seam_nodes(Model7(path))
    eP = sorted({v for e in _bedges(fP) for v in e})
    PE = np.array([X[v] for v in eP])
    pedge = sorted({eP[int(np.argmin(np.linalg.norm(PE - X[v], axis=1)))] for v in seam})
    ca = columns(f, (AVW, CX), 'SlidingElastic1Secondary', seam)
    cp = columns(f, (PVW, PEB), 'SlidingElastic1Primary', pedge)
    for name, cc in (('AVW / cervix', ca), ('PVW / PeB', cp)):
        L = Counter(len(c) for c in cc.values())
        th = [np.linalg.norm(f.nodes[c[-1]] - f.nodes[c[0]]) for c in cc.values()]
        print(f'{name}: {len(cc)} lumen-edge nodes; nodes per column {dict(sorted(L.items()))}; column length median '
              f'{np.median(th):.2f} mm (min {min(th):.2f}, max {max(th):.2f})')
    return f, seam, pedge, ca, cp


if __name__ == '__main__':
    main(*sys.argv[1:])
