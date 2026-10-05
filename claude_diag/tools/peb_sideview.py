"""Side view (y forward, z up) of the perineal body, the distal PVW/AVW and the LA near the midline, for several runs at
one time, with the body's best-fit rotation in each panel title.

usage: py -3.10 peb_sideview.py RUN [RUN ...] [--t 1.0] [--slab 6] --out PNG
  --slab: half-width in x (mm) of the midline slab drawn for the solids and the LA (default 6).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from febmodel import Feb  # noqa: E402
from xplt import Xplt  # noqa: E402
from peb_rotation import kabsch  # noqa: E402

COLS = (('_PickedSet347', 'tab:blue', 'AVW'), ('_PickedSet64', 'tab:red', 'PVW'), ('_PickedSet346', 'tab:orange', 'cervix'),
        ('_PickedSet66', 'k', 'perineal body'), ('LA_PCMPRM', 'tab:green', 'LA_PCMPRM'), ('LA_PCM', 'yellowgreen', 'LA_PCM'))


def main():
    args = sys.argv[1:]
    t_want, slab, out = 1.0, 6.0, None
    for flag in ('--t', '--slab', '--out'):
        if flag in args:
            i = args.index(flag)
            v = args[i + 1]
            args = args[:i] + args[i + 2:]
            if flag == '--t':
                t_want = float(v)
            elif flag == '--slab':
                slab = float(v)
            else:
                out = v
    runs = args
    fig, axs = plt.subplots(1, len(runs) + 1, figsize=(6.5 * (len(runs) + 1), 7), squeeze=False)
    axs = axs[0]
    for j, run in enumerate(runs):
        fe = Feb(os.path.join(RUNS_DIR, run, run + '.feb'))
        x = Xplt(os.path.join(RUNS_DIR, run, run + '.xplt'))
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        conv = [k for k, st in enumerate(x.states) if st[1] == 0]
        tt = np.array([x.states[k][0] for k in conv])
        k = conv[int(np.argmin(np.abs(tt - t_want)))]
        u = x.var(k, 'displacement')
        panels = [(axs[0], x.X, 'reference')] if j == 0 else []
        panels.append((axs[j + 1], x.X + u, f'{run}\nt = {x.states[k][0]:.3f}'))
        peb = np.array(sorted(idx[n] for n in fe.domain_nodes('_PickedSet66')))
        _, rv, _ = kabsch(x.X[peb], x.X[peb] + u[peb])
        for ax, xyz, title in panels:
            for dom, col, lab in COLS:
                if dom not in fe.elem_blocks:
                    continue
                ii = np.array(sorted(idx[n] for n in fe.domain_nodes(dom)))
                p = xyz[ii]
                sel = np.abs(x.X[ii][:, 0] - 1.4) < slab
                ax.scatter(p[sel, 1], p[sel, 2], s=2, c=col, label=lab)
            for ds, col in (('LA_sphincter_side_conn', 'm'), ('LA_sphincter_post_conn', 'c'), ('Parcus_conn', 'gold')):
                for a, b in fe.discsets.get(ds, []):
                    pa, pb = xyz[idx[a]], xyz[idx[b]]
                    ax.plot([pa[1], pb[1]], [pa[2], pb[2]], c=col, lw=0.4)
            ax.set_aspect('equal')
            ax.set_xlim(-40, 60)
            ax.set_ylim(-80, 0)
            ax.grid(alpha=0.3)
            ax.set_xlabel('y (forward)')
            ax.set_ylabel('z (up)')
            if title == 'reference':
                ax.set_title('reference')
                ax.legend(fontsize=7, loc='upper right', markerscale=4)
            else:
                ax.set_title(f'{title}; body rot_-x {-np.degrees(rv[0]):+.1f} deg', fontsize=9)
    plt.tight_layout()
    plt.savefig(out, dpi=75)
    print(out)


if __name__ == '__main__':
    main()
