"""Depth slices (lateral x AP, canal frame) through the tube's distal junction: every PVW / PeB / curve hex cut by the plane
depth = z, the curve hexes inside a PVW / PeB hex in cyan, the point-tied nodes within 0.6 mm of the plane in red.
usage: py -3.10 junction_slices.py RUN OUT.png [z ...] [--side 1|-1] [--at RUN T]"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load, frame, overlap, deformed  # noqa: E402
import re  # noqa: E402

EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
COL = {'_PickedSet64': ('#e8d23a', 'PVW'), '_PickedSet66': ('#e07070', 'PeB'), 'canal_wrap_hex_pvw': ('#b03cc0', 'curve (PVW side)'),
       'canal_wrap_hex_avw': ('#d94fd9', 'curve (AVW side)'), '_PickedSet347': ('#c98b5a', 'AVW'), 'over': ('#30d0ff', 'curve hex inside PVW / PeB')}


def cut(P, z):
    pts = []
    for a, b in EDGES:
        ha, hb = P[a, 2] - z, P[b, 2] - z
        if ha * hb < 0:
            pts.append(P[a] + (P[b] - P[a]) * ha / (ha - hb))
    if len(pts) < 3:
        return None
    Q = np.array(pts)[:, :2]
    c = Q.mean(0)
    return Q[np.argsort(np.arctan2(Q[:, 1] - c[1], Q[:, 0] - c[0]))]


def main(a):
    side = int(a[a.index('--side') + 1]) if '--side' in a else 1
    at = a[a.index('--at') + 1:a.index('--at') + 3] if '--at' in a else None
    rest = [v for i, v in enumerate(a) if v not in ('--side', '--at') and (i == 0 or a[i - 1] not in ('--side',))
            and not (at and v in at)]
    run, out, zs = rest[0], rest[1], [float(v) for v in rest[2:]] or [1.0, 2.5, 4.0, 6.0]
    t, X, E = load(run)
    c0, B = frame(X, E)
    over = overlap(X, E)
    title = run
    if at:
        X, tt = deformed(X, at[0], float(at[1])); title = f'{at[0]} at t {tt:.3f}'
    L = {k: (p - c0) @ B for k, p in X.items()}
    blk = re.search(r'<constraint name="(?:canal_wrap_embed|junction_ties)".*?</constraint>', t, re.S)
    tied = {int(m) for m in re.findall(r'<linear_constraint>\s*<node id="(\d+)" bc="x">', blk.group(0))} if blk else set()
    fig, axs = plt.subplots(1, len(zs), figsize=(4.2 * len(zs), 4.6))
    for ax, z in zip(np.atleast_1d(axs), zs):
        for k, (n, c) in E.items():
            if n not in COL:
                continue
            P = np.array([L[v] for v in c])
            if P[:, 0].mean() * side < 1 or P[:, 2].min() > z or P[:, 2].max() < z:
                continue
            Q = cut(P, z)
            if Q is None:
                continue
            key = 'over' if k in over else n
            ax.add_patch(Polygon(Q * [side, 1], closed=True, fc=COL[key][0], ec='k', lw=0.3, alpha=0.75 if n.startswith('canal') else 0.5))
        for v in tied:
            if abs(L[v][2] - z) < 0.6 and L[v][0] * side > 0:
                ax.plot(L[v][0] * side, L[v][1], 'o', color='red', ms=5, mec='k', mew=0.5)
        ax.set_title(f'depth {z:g} mm', fontsize=10)
        ax.set_xlim(2, 16); ax.set_ylim(-10, 4); ax.set_aspect('equal')
        ax.set_xlabel('lateral (mm, ' + ('right' if side > 0 else 'left, mirrored') + ')'); ax.grid(alpha=0.2)
    np.atleast_1d(axs)[0].set_ylabel('AP (mm, anterior +)')
    hs = [Polygon([[0, 0]], fc=c, ec='k', lw=0.3, label=l) for c, l in COL.values() if l != 'AVW' and 'AVW side' not in l]
    hs.append(plt.Line2D([], [], marker='o', ls='', color='red', mec='k', label='tied node (point ties now; facet ties in the rebuild)'))
    fig.legend(handles=hs, loc='lower center', ncol=5, fontsize=9, frameon=False)
    fig.suptitle(title + ': depth slices through the distal junction (depth = mm from the PM plane, hymen end)', fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == '__main__':
    main(sys.argv[1:])
