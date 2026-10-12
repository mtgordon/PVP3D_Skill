"""Element volume change and folding from a FEBio plot file, per hex8 domain group.

For every hex8 element of the domains whose name starts with PREFIX, the volume ratio V / V0 from the plotted
displacements (the hex split into 5 tets), at the last converged state (or the last one at or before --t):
the largest ratios (swelling / ballooning), how many exceed 1.5x and 2x, the smallest ratios, and the largest
ratio per 0.1 of time (when it starts). A NEGATIVE ratio means the hex is folded inside out somewhere (one of
its tets is inverted) even though FEBio's Gauss-point Jacobians were positive at convergence: the shape a user
sees as "self-penetration" or "explosion" in FEBio Studio.

usage: python element_volume_report.py model.feb model.xplt [model.log] --prefix NAME [--t 0.9]
(the .log is used to pick the converged states when the plot file has no status flags)"""
import sys

import numpy as np

from xplt_reader import Xplt, converged_state_indices

TETS = [(0, 1, 2, 5), (0, 2, 3, 7), (0, 5, 7, 4), (2, 7, 5, 6), (0, 2, 7, 5)]


def vols(X, conn):
    v = np.zeros(len(conn))
    for a, b, c, d in TETS:
        p = X[conn[:, a]]
        v += np.einsum('ij,ij->i', X[conn[:, b]] - p, np.cross(X[conn[:, c]] - p, X[conn[:, d]] - p)) / 6.0
    return v


def main(argv):
    tw, prefix, files = None, '', []
    while argv:
        a = argv.pop(0)
        if a == '--t':
            tw = float(argv.pop(0))
        elif a == '--prefix':
            prefix = argv.pop(0)
        else:
            files.append(a)
    xplt = next(f for f in files if f.endswith('.xplt'))
    log = next((f for f in files if f.endswith('.log')), None)
    x = Xplt(xplt)
    doms = [d for d in x.domains if str(d.get('name', '')).startswith(prefix) and np.asarray(d['conn']).shape[1] == 8]
    if not doms:
        sys.exit(f'no hex8 domain named {prefix}* in {xplt}: {[d.get("name") for d in x.domains]}')
    conn = np.vstack([np.asarray(d['conn']) for d in doms])
    eids = np.concatenate([np.asarray(d['eids']) for d in doms])
    V0 = vols(x.X, conn)
    idx = converged_state_indices(x, log)
    if tw is not None:
        idx = [i for i in idx if x.states[i][0] <= tw + 1e-9]
    hist = {}
    for i in idx[::max(1, len(idx) // 40)] + [idx[-1]]:
        t = x.states[i][0]
        r = vols(x.X + x.var(i, 'displacement'), conn) / V0
        b = round(np.floor(t * 10 + 1e-9) / 10, 1)
        hist[b] = max(hist.get(b, 0), r.max())
    print(f'{len(eids)} hex8 in {[d.get("name") for d in doms]}, t = {t:.3f}')
    print('largest V/V0 :', ', '.join(f'{eids[j]} {r[j]:.2f}' for j in np.argsort(-r)[:8]))
    print(f'> 2x: {(r > 2).sum()}   > 1.5x: {(r > 1.5).sum()}')
    print('smallest V/V0:', ', '.join(f'{eids[j]} {r[j]:.2f}' for j in np.argsort(r)[:8]), '(negative = folded)')
    print('largest by time:', ' '.join(f'{b:.1f}:{v:.1f}' for b, v in sorted(hist.items())))


if __name__ == '__main__':
    main(sys.argv[1:])
