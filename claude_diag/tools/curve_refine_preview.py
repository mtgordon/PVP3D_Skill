"""Picture of the proposed curve refinement (2026-10-15; the user: "could we remesh the curve so that the thickness stays the
same, but the other two dimensions are meshed more closely together? ... show me the picture first").
Curve hex local directions (L283): edge 0-1 = around the curve's U (~2 mm), 1-2 = through the wall (~0.6 mm, 4 layers = the
2.4 mm wall), 0-4 = along the canal (1.0-1.6 mm; the ring split of build_batch280 halved it within 18 mm of the PeB).
Proposal: every curve hex with its centre at rest z < ZCUT split in 2 AROUND the U (edges 0-1, 3-2, 4-5, 7-6 halved);
the through-wall layers are untouched. Left: one ring's cross-section (looking along the canal); right: the outer surface of
the distal zone on one side.
usage: py -3.10 curve_refine_preview.py RUN OUT.png [--side 1|-1] [--zcut -15] [--ringz -33]"""
import os
import sys
from collections import defaultdict, deque

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from feb_model import Feb  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
DOMS = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
RINGFACES = [(1, 2, 6, 5), (3, 0, 4, 7), (0, 1, 5, 4), (2, 3, 7, 6)]   # around / through-wall neighbours


def main(argv):
    opts = {'--side': 1.0, '--zcut': -15.0, '--ringz': -33.0}
    for k in list(opts):
        if k in argv:
            i = argv.index(k); opts[k] = float(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    run, out = argv[0], argv[1]
    f = Feb(os.path.join(RUNS, run, run + '.feb')); X = f.nodes
    H = {e: c for d in DOMS for e, c in f.elem_blocks[d][1].items()}
    cen = {e: np.mean([X[n] for n in c], 0) for e, c in H.items()}
    side = [e for e in H if np.sign(cen[e][0]) == opts['--side']]
    zone = [e for e in side if cen[e][2] < opts['--zcut']]
    # one ring near ringz: hexes joined through their around / through-wall faces
    byface = defaultdict(list)
    for e in side:
        for fc in RINGFACES:
            byface[frozenset(H[e][i] for i in fc)].append(e)
    seed = min(side, key=lambda e: abs(cen[e][2] - opts['--ringz']) + 0.01 * abs(cen[e][1]))
    ring, q = {seed}, deque([seed])
    while q:
        e = q.popleft()
        for fc in RINGFACES:
            for o in byface[frozenset(H[e][i] for i in fc)]:
                if o not in ring:
                    ring.add(o); q.append(o)
    ax_dir = np.mean([X[H[e][4]] - X[H[e][0]] for e in ring], 0); ax_dir /= np.linalg.norm(ax_dir)
    u = np.cross(ax_dir, [1.0, 0, 0]); u /= np.linalg.norm(u); v = np.cross(ax_dir, u)
    o = np.mean([cen[e] for e in ring], 0)
    P = lambda p: ((p - o) @ u, (p - o) @ v)  # noqa: E731

    fig = plt.figure(figsize=(15, 9))
    for k, split in enumerate((False, True)):
        a = fig.add_subplot(2, 2, 1 + k)
        for e in ring:
            c = H[e]; q4 = [X[c[i]] for i in (0, 1, 2, 3)]
            quads = [q4]
            if split:
                m01 = (q4[0] + q4[1]) / 2; m32 = (q4[3] + q4[2]) / 2
                quads = [[q4[0], m01, m32, q4[3]], [m01, q4[1], q4[2], m32]]
            for qq in quads:
                a.add_patch(Polygon([P(p) for p in qq], closed=True, fc='#cfe3f5' if split else '#f3dcc4', ec='k', lw=0.6))
        a.set_aspect('equal'); a.autoscale_view()
        a.set_title(('AFTER: ' if split else 'BEFORE: ') + f'one curve ring (rest z {o[2]:.0f}), looking along the canal\n'
                    + ('every element halved around the U; the 4 through-wall layers unchanged' if split else
                       f'{len(ring)} elements: 4 layers through the 2.4 mm wall x {len(ring) // 4} around the U'), fontsize=10)
        a.set_xlabel('mm'); a.set_ylabel('mm')

        b = fig.add_subplot(2, 2, 3 + k, projection='3d')
        polys = []
        bf = defaultdict(int)
        for e in side:
            for fc in RINGFACES[2:]:
                bf[frozenset(H[e][i] for i in fc)] += 1
        for e in zone:
            c = H[e]
            for fc, around in (((0, 1, 5, 4), ((0, 1), (4, 5))), ((2, 3, 7, 6), ((3, 2), (7, 6)))):
                if bf[frozenset(c[i] for i in fc)] != 1:
                    continue
                q4 = [X[c[i]] for i in fc]
                if split:
                    (a0, a1), (b0, b1) = around
                    m1 = (X[c[a0]] + X[c[a1]]) / 2; m2 = (X[c[b0]] + X[c[b1]]) / 2
                    if fc == (0, 1, 5, 4):
                        polys += [[X[c[0]], m1, m2, X[c[4]]], [m1, X[c[1]], X[c[5]], m2]]
                    else:
                        polys += [[X[c[2]], m1, m2, X[c[6]]], [m1, X[c[3]], X[c[7]], m2]]
                else:
                    polys.append(q4)
        b.add_collection3d(Poly3DCollection(polys, fc='#cfe3f5' if split else '#f3dcc4', ec='k', lw=0.3))
        allp = np.array([p for pp in polys for p in pp])
        for i, lab in enumerate('xyz'):
            getattr(b, f'set_{lab}lim')(allp[:, i].min(), allp[:, i].max())
        b.set_box_aspect(np.ptp(allp, 0))
        b.view_init(elev=20, azim=-35 if opts['--side'] > 0 else 215)
        b.set_title(('AFTER' if split else 'BEFORE') + f': curve surface below rest z {opts["--zcut"]:.0f} '
                    f'({"right" if opts["--side"] > 0 else "left"} side, seen from the side)\n'
                    + (f'{2 * len(zone)} elements' if split else f'{len(zone)} elements'), fontsize=10)
        b.set_xlabel('x'); b.set_ylabel('y'); b.set_zlabel('z')
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print('ring', len(ring), 'elements at rest z', round(o[2], 1), '; zone', len(zone), 'elements on this side ->', 2 * len(zone))


if __name__ == '__main__':
    main(sys.argv[1:])
