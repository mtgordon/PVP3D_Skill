"""Batch 452 (2026-10-23; the user: "build the next run with the AVW-Para ends at the AVW top moved from the curves onto the
nearest AVW nodes"). NOT IN SOURCE (spring ends): in AVW-Para-L_conn / AVW-Para-R_conn, every end within ZBAND mm of the AVW's
distal end (rest z) on the anterior side (rest y > YMIN) that is not an AVW node (a curve node or a spring-patch anchor on the
curve) moves to the nearest AVW (_PickedSet347) node not already taken by another moved end of the family. Spring material and
the other ends unchanged; old patch anchors stay (no spring on them).
usage: py -3.10 build_batch452.py SRC_RUN NEW_RUN [--zband 3] [--ymin 14]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

FAMS = ('AVW-Para-L_conn', 'AVW-Para-R_conn')


def main(src, new, zband, ymin):
    t, X, E = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    avw = sorted({v for n, c in E.values() if n == '_PickedSet347' for v in c})
    A = np.array([X[v] for v in avw])
    zmin = A[:, 2].min()
    report = []
    for fam in FAMS:
        m = re.search(r'(<DiscreteSet name="%s">)(.*?)(</DiscreteSet>)' % re.escape(fam), t, re.S)
        rows = re.findall(r'(\s*<delem>)(\d+),(\d+)(</delem>)', m.group(2))
        assert len(rows) == len(re.findall('<delem', m.group(2))), fam
        taken, moved, out = set(), [], []
        for o, p, q, c in rows:
            p, q = int(p), int(q)
            ends = []
            for v in (p, q):
                if v not in avw and X[v][2] < zmin + zband and X[v][1] > ymin:
                    d = np.linalg.norm(A - X[v], axis=1)
                    for i in np.argsort(d):
                        if avw[i] not in taken:
                            break
                    taken.add(avw[i]); moved.append((v, avw[i], round(float(d[i]), 2))); v = avw[i]
                ends.append(v)
            out.append(f'{o}{ends[0]},{ends[1]}{c}')
        body = ''.join(out) + re.search(r'\s*$', m.group(2)).group(0)
        t = t[:m.start(2)] + body + t[m.end(2):]
        report.append(f'{fam}: {len(moved)} ends moved (old -> AVW node, mm) {moved}')
    msg = ('NOT IN SOURCE (spring ends): AVW-Para ends at the AVW top (rest z < AVW distal end + %g, y > %g) moved from the curves '
           'onto the nearest AVW nodes: ' % (zband, ymin) + '; '.join(report))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch452.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    opt = lambda n, d: (a[a.index(n) + 1] if n in a else d)  # noqa: E731
    main(a[0], a[1], float(opt('--zband', 3)), float(opt('--ymin', 14)))
