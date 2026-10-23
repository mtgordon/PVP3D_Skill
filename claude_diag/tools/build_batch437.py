"""Batch 437 (2026-10-23; L437 failed at its first step with negative jacobians even at dt 3e-4: after the distal curve splits,
22 curve_all nodes sit 0-0.49 mm BEHIND another curve_all face at rest (L423b: 0), around the strap rim corner (split edge
midpoints on the chord inside the curved surface) -> the self-contact kicks them at the first step).
NOT IN SOURCE (our contact): every face of SURF that has a node of SURF (not its own) behind it within SR at rest (projecting
inside the face) is removed from SURF; the rest stays. As build_batch340's rest check.
usage: py -3.10 build_batch437.py SRC_RUN NEW_RUN [--surf curve_all] [--sr 0.5]"""
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def behind(X, f, p, sr):
    P = [X[q] for q in f]; A = P[0]
    n = np.cross(P[1] - P[0], P[2] - P[0]); n /= np.linalg.norm(n)
    dd = (p - A) @ n
    if not (-sr < dd < 0.0):
        return False
    q = p - dd * n
    for a, b, c in [(P[0], P[1], P[2])] + ([(P[0], P[2], P[3])] if len(P) == 4 else []):
        v0, v1, v2 = b - a, c - a, q - a
        d00, d01, d11, d20, d21 = v0 @ v0, v0 @ v1, v1 @ v1, v2 @ v0, v2 @ v1
        den = d00 * d11 - d01 * d01
        l1 = (d11 * d20 - d01 * d21) / den; l2 = (d00 * d21 - d01 * d20) / den
        if l1 >= -0.01 and l2 >= -0.01 and l1 + l2 <= 1.01:
            return True
    return False


def main(src, new, surf, sr):
    t, X, E = load(src)
    X = {k: np.asarray(v, float) for k, v in dict(X).items()}
    m = re.search(r'(<Surface name="%s">\n)(.*?)(\t\t</Surface>)' % surf, t, re.S)
    rows = re.findall(r'<(quad4|tri3) id="\d+">([\d,]+)</', m.group(2))
    faces = [(k, list(map(int, f.split(',')))) for k, f in rows]
    nodes = sorted({v for _, f in faces for v in f})
    C = np.array([np.mean([X[v] for v in f], 0) for _, f in faces]); tr = cKDTree(C)
    drop = set()
    for v in nodes:
        for j in tr.query_ball_point(X[v], sr + 3):
            if v not in faces[j][1] and behind(X, faces[j][1], X[v], sr):
                drop.add(j)
    keep = [faces[j] for j in range(len(faces)) if j not in drop]
    body = ''.join(f'\t\t\t<{k} id="{i}">' + ','.join(map(str, f)) + f'</{k}>\n' for i, (k, f) in enumerate(keep, 1))
    t = t[:m.start(2)] + body + t[m.end(2):]
    zs = [np.mean([X[v][2] for v in faces[j][1]]) for j in drop]
    msg = (f'NOT IN SOURCE (contact surface): {len(drop)} of {len(faces)} {surf} faces with a {surf} node behind them within {sr:g} mm '
           f'at rest removed (rest z {min(zs):.1f}..{max(zs):.1f})' if drop else f'{surf}: no face with a node behind it')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch437.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], a[a.index('--surf') + 1] if '--surf' in a else 'curve_all', float(a[a.index('--sr') + 1]) if '--sr' in a else 0.5)
