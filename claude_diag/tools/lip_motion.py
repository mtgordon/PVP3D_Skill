"""How do the seam-lip nodes move late in a run? (2026-10-08, the user: should the lips be stiff in the direction of the force
and soft in the others?) The lips (tools/seam_extend.py: hexes AVW_ext / cervix_ext / PVW_ext, each built as a wall side face
'lo' + its outward copy 'hi') get a local frame per element: D = outward (mean hi - lo), A = along the canal (the side face's
edge most perpendicular to D and the lip normal... taken as the face edge closest to the canal axis), N = D x A (the lip normal,
through the lip). Each lip node averages the frames of its elements. Over the last N plot states, every step of every lip node
is split into its N / D / A parts; reported: share of the squared step in each direction for the largest steps, and how the
step size grows with the ring (distance from the seam: 0 = wall node, 1..K = lip rings).
usage: py -3.10 lip_motion.py RUN [N]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
LIPS = ('AVW_ext', 'cervix_ext', 'PVW_ext')

if __name__ == '__main__':
    run = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    X0 = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(',')])
          for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t)}
    lipnodes = {int(m.group(1)) for m in re.finditer(r'<node id="(\d+)">',
                re.search(r'<Nodes name="seam_ext_nodes">(.*?)</Nodes>', t, re.S).group(1))}
    # canal axis: principal axis of the AVW nodes
    avw = re.search(r'<Elements type="hex8" name="_PickedSet347">(.*?)</Elements>', t, re.S).group(1)
    P = np.array([X0[int(v)] for e in re.findall(r'>([\d,]+)</elem>', avw) for v in e.split(',')])
    ax = np.linalg.svd(P - P.mean(0), full_matrices=False)[2][0]
    frames, owner = {}, {}
    for dn in LIPS:
        blk = re.search(r'<Elements type="hex8" name="%s">(.*?)</Elements>' % dn, t, re.S).group(1)
        for e in re.findall(r'>([\d,]+)</elem>', blk):
            h = [int(v) for v in e.split(',')]
            lo, hi = np.array([X0[v] for v in h[:4]]), np.array([X0[v] for v in h[4:]])
            D = (hi - lo).mean(0); D /= np.linalg.norm(D)
            edges = [lo[(k + 1) % 4] - lo[k] for k in range(4)]
            A = max(edges, key=lambda v: abs(v @ ax) / np.linalg.norm(v)); A = A - (A @ D) * D; A /= np.linalg.norm(A)
            N = np.cross(D, A)
            for v in h:
                frames.setdefault(v, []).append((D, A, N)); owner[v] = dn
    F = {}
    for v, fs in frames.items():
        D = np.mean([f[0] for f in fs], 0); D /= np.linalg.norm(D)
        A = np.mean([f[1] for f in fs], 0); A -= (A @ D) * D; A /= np.linalg.norm(A)
        F[v] = (D, A, np.cross(D, A))
    # ring: 0 = wall node, else by distance to the nearest wall node of the lip
    wallv = [v for v in F if v not in lipnodes]
    W = np.array([X0[v] for v in wallv])
    dist = {v: (0.0 if v not in lipnodes else float(np.min(np.linalg.norm(W - X0[v], axis=1)))) for v in F}
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    idx = {int(k): i for i, k in enumerate(x.node_ids)}
    ids = list(range(max(0, len(x.states) - n - 1), len(x.states)))
    vs = sorted(F)
    prev = None
    rows = []
    for i in ids:
        u = np.asarray(x.var(i, 'displacement'), dtype=float)[[idx[v] for v in vs]]
        if prev is not None:
            du = u - prev
            for j, v in enumerate(vs):
                D, A, N = F[v]
                rows.append((x.states[i][0], v, np.linalg.norm(du[j]), du[j] @ N, du[j] @ D, du[j] @ A))
        prev = u
    R = np.array([r[2:] for r in rows])
    big = R[:, 0] >= np.percentile(R[:, 0], 99)
    def share(M):
        s = (M[:, 1:] ** 2).sum(0)
        return s / s.sum()
    print(f'{run}: {len(vs)} lip-element nodes ({len(lipnodes)} lip, {len(vs) - len(lipnodes)} wall), last {len(ids) - 1} steps')
    for lab, M in (('all steps', R), ('largest 1 % of steps', R[big])):
        sN, sD, sA = share(M)
        print(f'  {lab}: share of squared motion  normal (through the lip) {sN:.0%}  outward (in-plane) {sD:.0%}  along the canal {sA:.0%}')
    # growth with the distance from the seam (flap = grows with distance)
    d = np.array([dist[r[1]] for r in rows])
    for lo_, hi_ in ((-0.1, 0.1), (0.1, 2.0), (2.0, 3.6), (3.6, 6.0)):
        m = (d > lo_) & (d <= hi_)
        if m.any():
            print(f'  distance from the seam {max(lo_, 0):.1f}-{hi_:.1f} mm: median step {np.median(R[m, 0]):.3f} mm, p99 '
                  f'{np.percentile(R[m, 0], 99):.3f}, median |normal part| {np.median(np.abs(R[m, 1])):.3f}')
    # sign alternation of the normal part for the 5 most active lip nodes (back-and-forth = flapping)
    act = {}
    for r in rows:
        act.setdefault(r[1], []).append(r[3])
    top = sorted(act, key=lambda v: -np.sum(np.abs(act[v])))[:5]
    for v in top:
        s = np.sign(act[v]); flips = int(np.sum(s[1:] * s[:-1] < 0))
        print(f'  node {v} ({owner[v]}, {dist[v]:.1f} mm out): normal steps sum {np.sum(act[v]):+.2f}, |sum| of steps '
              f'{np.sum(np.abs(act[v])):.2f} mm, sign flips {flips} / {len(s) - 1}')
