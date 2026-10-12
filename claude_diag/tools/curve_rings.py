"""The tube's curves as structured rings (2026-10-09). Each curve pair (canal_wrap_hex_avw + _pvw, one U per lateral side) is a
stack of 5 x 7 node rings (5 through the thickness, row 0 = the U's inner surface; 7 round the U, column 0 on the AVW / cervix
side, column 6 on the PVW side, column 3 = where the AVW half meets the PVW half). Walked from the apical end cap (the same
cap apex_contour.py starts from) down the canal through the hexes. usage (module): rings(t) -> {sgn: [grid, ...]} (apex first)"""
import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apex_contour import sweep  # noqa: E402

CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')


def rings(t):
    sides, c0, B, X = sweep(t)
    mesh = ET.fromstring(t).find('Mesh')
    H = [[int(v) for v in el.text.split(',')] for e in mesh.findall('Elements') if e.get('name') in CURVES for el in e]
    byface = defaultdict(list)
    for h in H:
        byface[frozenset(h[:4])].append(h); byface[frozenset(h[4:])].append(h)
    out = {}
    for (Gr, R, msg), sgn in zip(sides, (1, -1)):
        if (X[Gr[0, 0]] - c0) @ B[:, 1] < (X[Gr[0, 6]] - c0) @ B[:, 1]:
            Gr = Gr[:, ::-1]                                # column 0 = the anterior (AVW / cervix) side
        grids, used = [Gr], set()
        while True:
            G = grids[-1]; N = -np.ones_like(G)
            for a in range(4):
                for b in range(6):
                    q = [G[a, b], G[a, b + 1], G[a + 1, b + 1], G[a + 1, b]]
                    hs = [h for h in byface[frozenset(q)] if tuple(h) not in used]
                    if not hs:
                        continue
                    h = hs[0]; used.add(tuple(h))
                    for v in q:
                        i = h.index(v); o = h[(i + 4) % 8]
                        N[np.where(G == v)] = o
            if (N < 0).any():
                break
            grids.append(N)
        out[sgn] = (grids, c0, B, X)
    return out


if __name__ == '__main__':
    from variants7 import JOBS
    run = sys.argv[1]
    t = open(os.path.join(JOBS, 'claude_diag', 'runs', run, run + '.feb'), encoding='utf-8').read()
    for sgn, (grids, c0, B, X) in rings(t).items():
        L = lambda v: (X[v] - c0) @ B  # noqa: E731
        print(sgn, len(grids), 'rings; depth', round(L(grids[0][2, 3])[2], 1), '->', round(L(grids[-1][2, 3])[2], 1))
        for g in (grids[0], grids[len(grids) // 2]):
            for a in (0, 4):
                print('  row', a, ' '.join('(%.1f,%.1f,%.1f)' % tuple(L(v)) for v in g[a]))
