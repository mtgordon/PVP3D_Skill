"""Picture of the proposed spring-end patches (2026-10-15; the user: "I don't like how the springs ... attaching to the material
only pull at one node. Is there a way to distribute the force to the line instead"; "show me the patches first").
For every discrete spring end on a curve node at rest z < ZCUT: the patch = the curve's boundary faces around that node
(one ring); weights = each face's area / 4 to each of its corners, normalised (a bilinear 'tent': on a regular grid the
centre node 1/4, the 4 edge neighbours 1/8, the 4 corner neighbours 1/16). The spring would end on a new anchor node tied to
the weighted average of the patch (linear constraints = a distributing coupling).
usage: py -3.10 spring_patch_preview.py RUN OUT.png [--side 1|-1] [--zcut -15]"""
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from feb_model import Feb  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def curve_faces(f):
    H = [c for name, (et, d) in f.elem_blocks.items() if name.startswith('canal_wrap_hex') and name != 'canal_wrap_hex_avw_cervix'
         for c in d.values()]
    cnt = Counter(frozenset(h[i] for i in fc) for h in H for fc in HF)
    return [[h[i] for i in fc] for h in H for fc in HF if cnt[frozenset(h[i] for i in fc)] == 1]


def _normal(X, q):
    P = np.array([X[v] for v in q])
    n = np.cross(P[2] - P[0], P[3] - P[1])
    return n / np.linalg.norm(n)


def _axial(f):
    """node -> neighbours along the canal (curve hex edges 0-4, 1-5, 2-6, 3-7)"""
    if not hasattr(f, '_ax'):
        ax = defaultdict(set)
        for name, (et, d) in f.elem_blocks.items():
            if name.startswith('canal_wrap_hex') and name != 'canal_wrap_hex_avw_cervix':
                for c in d.values():
                    for i in range(4):
                        ax[c[i]].add(c[i + 4]); ax[c[i + 4]].add(c[i])
        f._ax = ax
    return f._ax


def patches(f, zcut=-15.0, corner=False, corner_deg=60.0, keep_deg=45.0):
    # corner='edge' (2026-10-15, the user: the one-face tip patches pull off-centre, E29216): at a tip corner the force goes
    # along the corner edge (crease) instead: the node 1/2, its neighbours along the crease 1/4 each (3/4 + 1/4 at an end)
    """{end node: [(node, weight)]}, {end node: (spring set, other end)}
    corner=True (2026-10-15, the user saw the tip patches shear the corner): where the ring's faces meet at more than
    corner_deg (the curve's tip corners), keep only the faces whose outward normal is within keep_deg of the spring direction
    (node -> other end; mean over the node's springs), or the single best one if none is."""
    faces = curve_faces(f)
    X = f.nodes
    byn = defaultdict(list)
    for q in faces:
        for n in q:
            byn[n].append(q)
    ends = {}
    for name, pairs in f.discsets.items():
        for a, b in pairs:
            for n, o in ((a, b), (b, a)):
                if n in byn and X[n][2] < zcut:
                    ends.setdefault(n, []).append((name, o))
    out = {}
    for n in ends:
        ring = byn[n]
        # corner='node' (2026-10-15 night): the tip corner ends are left out of P (their springs stay on the node)
        if corner:
            N = [_normal(X, q) for q in ring]
            ang = max(np.degrees(np.arccos(np.clip(a @ b, -1, 1))) for a in N for b in N)
            if ang > corner_deg and corner == 'axial':   # 2026-10-16 (L288f): up the curve, no overlap with other tip ends
                ax = _axial(f)
                up1 = [v for v in ax.get(n, ()) if X[v][2] > X[n][2]] or list(ax.get(n, ()))
                if up1:
                    u1 = max(up1, key=lambda v: X[v][2])
                    up2 = [v for v in ax.get(u1, ()) if v != n and X[v][2] > X[u1][2]]
                    out[n] = [(n, 0.5), (u1, 1 / 3)] + ([(max(up2, key=lambda v: X[v][2]), 1 / 6)] if up2 else [])
                    if not up2:
                        out[n] = [(n, 0.5), (u1, 0.5)]
                else:
                    out[n] = [(n, 1.0)]
                continue
            if ang > corner_deg and corner == 'node':   # tip corner ends keep their single node (no patch)
                continue
            if ang > corner_deg and corner == 'edge':
                nb = []
                for i in range(len(ring)):
                    for j in range(i + 1, len(ring)):
                        sh = set(ring[i]) & set(ring[j])
                        if n in sh and len(sh) == 2 and np.degrees(np.arccos(np.clip(N[i] @ N[j], -1, 1))) > corner_deg:
                            nb.append((sh - {n}).pop())
                nb = sorted(set(nb))
                out[n] = ([(n, 0.5)] + [(v, 0.25) for v in nb[:2]]) if len(nb) >= 2 else                          ([(n, 0.75), (nb[0], 0.25)] if nb else [(n, 1.0)])
                continue
            if ang > corner_deg:
                d = np.mean([(X[o] - X[n]) / np.linalg.norm(X[o] - X[n]) for _, o in ends[n]], 0); d /= np.linalg.norm(d)
                dots = [nv @ d for nv in N]
                keep = [q for q, c in zip(ring, dots) if c >= np.cos(np.radians(keep_deg))]
                ring = keep or [ring[int(np.argmax(dots))]]
        w = defaultdict(float)
        for q in ring:
            P = np.array([X[v] for v in q])
            area = 0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[1]))
            for v in q:
                w[v] += area / 4
        s = sum(w.values())
        out[n] = [(v, x / s) for v, x in w.items()]
    return out, ends, faces


def main(argv):
    side, zcut = 1.0, -15.0
    if '--side' in argv:
        i = argv.index('--side'); side = float(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    if '--zcut' in argv:
        i = argv.index('--zcut'); zcut = float(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    run, out = argv[0], argv[1]
    f = Feb(os.path.join(RUNS, run, run + '.feb')); X = f.nodes
    P, ends, faces = patches(f, zcut)
    mine = {n: p for n, p in P.items() if np.sign(X[n][0]) == side}
    zone = [q for q in faces if np.sign(np.mean([X[v][0] for v in q])) == side and np.mean([X[v][2] for v in q]) < zcut + 3]
    fams = sorted({nm for n in mine for nm, _ in ends[n]})
    col = dict(zip(fams, plt.cm.tab10(np.linspace(0, 1, 10))[:len(fams)]))
    fig = plt.figure(figsize=(16, 9))
    for k, (elev, azim) in enumerate(((15, -40 if side > 0 else 220), (15, 40 if side > 0 else 140))):
        ax = fig.add_subplot(1, 2, k + 1, projection='3d')
        ax.add_collection3d(Poly3DCollection([[X[v] for v in q] for q in zone], fc=(0.93, 0.93, 0.93, 0.5), ec=(0.5, 0.5, 0.5, 0.4), lw=0.2))
        for n, p in mine.items():
            fam = ends[n][0][0]
            pf = [q for q in zone if n in q]
            ax.add_collection3d(Poly3DCollection([[X[v] for v in q] for q in pf], fc=(*col[fam][:3], 0.55), ec='k', lw=0.4))
            ax.scatter(*X[n], color='k', s=10)
            for nm, o in ends[n]:
                d = X[o] - X[n]; d = d / np.linalg.norm(d) * 4
                ax.plot(*zip(X[n], X[n] + d), color=col[nm], lw=1.2)
        allp = np.array([X[v] for q in zone for v in q])
        for i, lab in enumerate('xyz'):
            getattr(ax, f'set_{lab}lim')(allp[:, i].min(), allp[:, i].max())
        ax.set_box_aspect(np.ptp(allp, 0)); ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    hs = [plt.Line2D([], [], color=col[fm], lw=6) for fm in fams]
    fig.legend(hs, [f'{fm} ({sum(1 for n in mine if ends[n][0][0] == fm)} ends)' for fm in fams], loc='lower center', ncol=4, fontsize=9)
    nw = [len(p) for p in mine.values()]
    fig.suptitle(f'{run}: {"right" if side > 0 else "left"} curve below rest z {zcut:.0f}: each spring end (black dot, short line = '
                 f'spring direction) would pull on its coloured patch\n({len(mine)} ends; {min(nw)}-{max(nw)} nodes per patch, '
                 f'tent weights: centre ~1/4, edge neighbours ~1/8, corner neighbours ~1/16) instead of on the one node', fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    fig.savefig(out, dpi=120)
    print(len(mine), 'ends on this side; patch sizes', Counter(nw), '; families', {fm: sum(1 for n in mine if ends[n][0][0] == fm) for fm in fams})


if __name__ == '__main__':
    main(sys.argv[1:])
