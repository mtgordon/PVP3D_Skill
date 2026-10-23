"""How far the canal seam's edges slide past each other (2026-10-07, to size an unloaded wall extension at the seam).
The seam: SlidingElastic1's secondary (AVW + cervix) lateral edge nodes lie within 0.75 mm of the primary's (PVW + PeB)
boundary edge at rest (build_batch90.seam_nodes). Per plot state, for each secondary seam node: its OVERHANG past the
primary's current boundary edge, measured in the plane of the primary's edge face, perpendicular to the edge (+ = beyond
the primary's last facet, - = still over it). And the reverse: each primary boundary node within 0.75 mm of the secondary's
edge at rest against the secondary's current edge. Prints max / p95 / count beyond 0 / 1 / 2 / 3 / 5 mm.
usage: py -3.10 seam_slide.py RUN [EVERY]   (EVERY: use every n-th plot state, default 5; the last state always)"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402
from variants7 import Model7, JOBS  # noqa: E402
from build_batch79 import _bedges, _dseg  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')


def edge_faces(faces):
    """boundary edge -> the face it belongs to"""
    own = {}
    for f in faces:
        for k in range(len(f)):
            e = tuple(sorted((f[k], f[(k + 1) % len(f)])))
            own.setdefault(e, []).append(f)
    return {e: fs[0] for e, fs in own.items() if len(fs) == 1}


def overhang(P, edges, X):
    """signed in-plane distance of point P past the nearest boundary edge (edges: {(a, b): face})"""
    best = None
    for (a, b), f in edges.items():
        A, B = X[a], X[b]
        t = np.clip((P - A) @ (B - A) / ((B - A) @ (B - A)), 0, 1)
        Q = A + t * (B - A)
        d = np.linalg.norm(P - Q)
        if best is None or d < best[0]:
            best = (d, a, b, f, Q)
    _, a, b, f, Q = best
    e = X[b] - X[a]
    e /= np.linalg.norm(e)
    c = np.mean([X[v] for v in f], axis=0) - Q
    inward = c - (c @ e) * e
    inward /= np.linalg.norm(inward)
    return -(P - Q) @ inward


if __name__ == '__main__':
    run = sys.argv[1]
    every = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    m = Model7(os.path.join(RUNS, run, run + '.feb'))
    surf = {s.get('name'): s for s in m.mesh.findall('Surface')}
    fl = lambda s: [tuple(int(v) for v in f.text.split(',')) for f in s]  # noqa: E731
    fS, fP = fl(surf['SlidingElastic1Secondary']), fl(surf['SlidingElastic1Primary'])
    eS, eP = edge_faces(fS), edge_faces(fP)
    X0 = m.nodes()
    sec = sorted({v for e in eS for v in e if min(_dseg(X0[v], X0[a], X0[b]) for a, b in eP) < 0.75})
    pri = sorted({v for e in eP for v in e if min(_dseg(X0[v], X0[a], X0[b]) for a, b in eS) < 0.75})
    # only the nearby edges matter for each node: pre-select the 40 nearest edges at rest
    def near(v, edges):
        ks = list(edges)
        mid = np.array([(X0[a] + X0[b]) / 2 for a, b in ks])
        return {ks[i]: edges[ks[i]] for i in np.argsort(np.linalg.norm(mid - X0[v], axis=1))[:40]}
    nS = {v: near(v, eP) for v in sec}
    nP = {v: near(v, eS) for v in pri}
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    X0a = np.array([X0[int(n)] for n in x.node_ids])
    sts = list(range(0, len(x.states), every))
    if sts[-1] != len(x.states) - 1:
        sts.append(len(x.states) - 1)
    print(f'{run}: {len(sec)} secondary (AVW / cervix) seam nodes, {len(pri)} primary (PVW / PeB) seam nodes')
    print('    t   | sec past PVW edge: max  p95  n>0 n>1 n>2 n>3 n>5 | pri past AVW edge: max  p95  n>0 n>1 n>2 n>3 n>5')
    worst = {}
    for i in sts:
        u = np.asarray(x.var(i, 'displacement'), dtype=float)
        Xa = X0a + u
        X = {n: Xa[idx[n]] for n in set(sec) | set(pri) | {v for d in list(nS.values()) + list(nP.values()) for e in d for v in e}
             | {v for d in list(nS.values()) + list(nP.values()) for f in d.values() for v in f}}
        oS = np.array([overhang(X[v], nS[v], X) for v in sec])
        oP = np.array([overhang(X[v], nP[v], X) for v in pri])
        for v, o in zip(sec + pri, list(oS) + list(oP)):
            worst[v] = max(worst.get(v, -1e9), o)
        cnt = lambda o: ' '.join(f'{int((o > k).sum()):3d}' for k in (0, 1, 2, 3, 5))  # noqa: E731
        print(f' {x.states[i][0]:.4f} | {oS.max():20.2f} {np.percentile(oS, 95):4.1f} {cnt(oS)} | '
              f'{oP.max():20.2f} {np.percentile(oP, 95):4.1f} {cnt(oP)}')
    top = sorted(worst.items(), key=lambda kv: -kv[1])[:8]
    print('largest overhang over the run (node: mm, s = secondary, p = primary):',
          ', '.join(f'{v}{"s" if v in sec else "p"}: {o:.1f}' for v, o in top))
