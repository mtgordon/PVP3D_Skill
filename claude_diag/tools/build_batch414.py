"""Batch 414 (2026-10-22; the user: "I think the interaction between the outer curves has too many faces. It only needs to have
the outer faces and the front opening faces I believe"). NOT IN SOURCE (the curves and their self-contact are ours).
The curve self-contact surface curve_all (primary = secondary of curve_all_self) is cut to the curves' OUTER skin and their
DISTAL end (opening) faces. Each curve_all quad is classed by its hex face: 0123 / 4567 = ring ends (distal if rest z < -30,
else proximal); 0154 / 2376 = the U's skins (outer = the one farther in x-y from the canal centre at that level, the centre =
mean x-y of all curve nodes within 3 mm in z; inner = the other); 1265 / 3047 = the U's edges. Kept: outer + end_distal.
usage: py -3.10 build_batch414.py SRC_RUN NEW_RUN [--keep outer,end_distal]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

FT = {(0, 1, 2, 3): 'end', (4, 5, 6, 7): 'end', (0, 1, 5, 4): 'skinA', (2, 3, 7, 6): 'skinB',
      (1, 2, 6, 5): 'edge', (3, 0, 4, 7): 'edge'}


def classify(t, X, E, surface='curve_all'):
    cur = {k: c for k, (n, c) in E.items() if n.startswith('canal_wrap') and len(c) == 8}
    f2 = defaultdict(list)
    for k, c in cur.items():
        for f, ty in FT.items():
            f2[frozenset(c[i] for i in f)].append((k, ty))
    allc = np.array([X[v] for c in cur.values() for v in c])
    m = re.search(r'<Surface name="%s">(.*?)</Surface>' % surface, t, re.S)
    out = []
    for kind, fs in re.findall(r'<(quad4|tri3) id="\d+">([\d,]+)</', m.group(1)):
        f = [int(v) for v in fs.split(',')]
        hits = f2.get(frozenset(f))
        if not hits:
            out.append(('unknown', kind, f)); continue
        k, ty = hits[0]; c = cur[k]
        P = np.mean([X[v] for v in f], 0)
        if ty == 'end':
            out.append(('end_distal' if P[2] < -30 else 'end_prox', kind, f)); continue
        if ty in ('skinA', 'skinB'):
            other = (2, 3, 7, 6) if ty == 'skinA' else (0, 1, 5, 4)
            Q = np.mean([X[c[i]] for i in other], 0)
            cc = allc[np.abs(allc[:, 2] - P[2]) < 3][:, :2].mean(0)
            out.append(('outer' if np.linalg.norm(P[:2] - cc) > np.linalg.norm(Q[:2] - cc) else 'inner', kind, f)); continue
        out.append(('edge', kind, f))
    return out


def main(src, new, keep=('outer', 'end_distal')):
    t, X, E = load(src)
    X = {k: np.asarray(p, float) for k, p in dict(X).items()}
    cls = classify(t, X, E)
    n0 = Counter(c for c, _, _ in cls)
    kept = [(kind, f) for c, kind, f in cls if c in keep]
    body = ''.join(f'\t\t\t<{kind} id="{i}">' + ','.join(map(str, f)) + f'</{kind}>\n' for i, (kind, f) in enumerate(kept, 1))
    t = re.sub(r'(<Surface name="curve_all">\n).*?(\t\t</Surface>)', lambda m: m.group(1) + body + m.group(2), t, count=1, flags=re.S)
    msg = (f'NOT IN SOURCE: curve self-contact surface curve_all cut to {"+".join(keep)} faces (the user): {len(kept)} of '
           f'{len(cls)} facets kept; classes before {dict(n0)}')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch414.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    kp = tuple(a[a.index('--keep') + 1].split(',')) if '--keep' in a else ('outer', 'end_distal')
    main(a[0], a[1], kp)
