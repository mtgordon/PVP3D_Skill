"""Batch 361 (2026-10-20 night; the user's go (1) of NEXT_SESSION_PROMPT_2026-10-20b: "try including more nodes in that contact
set"). curve_all (the curves' self-contact surface, curve_all_self) holds 1632 of the 2696 curve boundary faces; the other 1064
all touch the wall attachment: 134 with every node tied (canal_wrap_ends / junction_ties dependent nodes, or nodes shared with a
wall / PeB element) and 930 partly tied. In L350b at t 0.70 16 of 111 penetrating curve nodes were on such faces.
NOT IN SOURCE (our self-contact): the PARTLY-tied curve boundary faces with every node at rest z < ZCUT (the collapse zone, rest
z -40 to -15) are added to curve_all, wound outward from their hex; a face is left out if any of its nodes lies within SR (the
search radius, 0.5 mm) of a curve_all / added face that shares no node with that node's own hexes (skill gotcha 30: the contact
must not grab the part's own opposite face at rest). Fully tied faces stay out (they cannot move apart from the wall).
usage: py -3.10 build_batch361.py SRC_RUN NEW_RUN [--zcut -12] [--sr 0.5] [--dry]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from build_batch162 import write  # noqa: E402
from feb_model import pt_tri_dist  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
FI = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def elements(t):
    out = {}
    for m in re.finditer(r'\t\t<Elements type="(\w+)" name="([^"]+)">\n(.*?)\t\t</Elements>\n', t, re.S):
        out[m.group(2)] = (m.group(1), [[int(v) for v in c.split(',')] for c in re.findall(r'>([0-9,]+)</elem>', m.group(3))])
    return out


def tied_nodes(t):
    dep = set()
    for nm in ('canal_wrap_ends', 'junction_ties'):
        blk = re.search(rf'<constraint name="{nm}".*?</constraint>', t, re.S).group(0)
        for lc in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', blk, re.S):
            dep.add(int(re.search(r'<node id="(\d+)" bc="\w">1</node>', lc).group(1)))
    return dep


def main(src, new, zcut, sr, dry):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    E = elements(t)
    curve = [d for d in E if d.startswith('canal_wrap_hex')]
    H = [h for d in curve for h in E[d][1]]
    cn = set(v for h in H for v in h)
    other = set(v for d, (ty, el) in E.items() if d not in curve for e in el for v in e)
    tied = tied_nodes(t) & cn   # constraint-tied only (the 10-19 count: 134 full / 930 partly)
    cnt = Counter(tuple(sorted(h[i] for i in f)) for h in H for f in FI)
    bfaces = []   # (quad, hex index)
    for k, h in enumerate(H):
        for f in FI:
            q = [h[i] for i in f]
            if cnt[tuple(sorted(q))] == 1:
                bfaces.append((q, k))
    ms = re.search(r'\t\t<Surface name="curve_all">\n(.*?)\t\t</Surface>\n', t, re.S)
    have = {tuple(sorted(int(v) for v in c.split(','))) for c in re.findall(r'>([0-9,]+)</quad4>', ms.group(1))}
    nfree = sum(1 for q, _ in bfaces if not any(v in tied for v in q))
    full = [(q, k) for q, k in bfaces if all(v in tied for v in q)]
    part = [(q, k) for q, k in bfaces if tuple(sorted(q)) not in have and not all(v in tied for v in q) and any(v in tied for v in q)]
    print(f'curve boundary faces {len(bfaces)}: in curve_all {sum(1 for q, _ in bfaces if tuple(sorted(q)) in have)}; free {nfree}; '
          f'fully tied {len(full)}; partly tied, not in curve_all {len(part)}')
    cand = [(q, k) for q, k in part if all(X[v][2] < zcut for v in q)]
    print(f'partly tied with every node at rest z < {zcut:g}: {len(cand)}')
    # winding: outward from the hex centre
    def wind(q, k):
        c = np.mean([X[v] for v in H[k]], 0); P = [X[v] for v in q]
        n = np.cross(P[2] - P[0], P[3] - P[1])
        return q if n @ (np.mean(P, 0) - c) > 0 else q[::-1]
    cand = [(wind(q, k), k) for q, k in cand]
    # rest check: every node of a candidate face vs every curve_all / candidate face sharing no node with that node's hexes
    allq = [[int(v) for v in c.split(',')] for c in re.findall(r'>([0-9,]+)</quad4>', ms.group(1))] + [q for q, _ in cand]
    tris = [(q[0], q[1], q[2]) for q in allq] + [(q[0], q[2], q[3]) for q in allq]
    tc = cKDTree(np.array([np.mean([X[v] for v in tr], 0) for tr in tris]))
    nhex = defaultdict(set)
    for h in H:
        for v in h:
            nhex[v] |= set(h)
    keep, drop = [], []
    for q, k in cand:
        dmin = 9e9
        for v in q:
            for j in tc.query_ball_point(X[v], sr + 6.0):
                tr = tris[j]
                if set(tr) & nhex[v]:
                    continue
                dmin = min(dmin, pt_tri_dist(X[v], *(X[u] for u in tr)))
        (drop if dmin < sr else keep).append((q, dmin))
    print(f'rest check (sr {sr:g}): keep {len(keep)}, drop {len(drop)}; nearest non-neighbour face of a kept face '
          f'{min(d for _, d in keep):.2f} mm' if keep else 'none kept')
    zs = [X[v][2] for q, _ in keep for v in q]; xs = [X[v][0] for q, _ in keep for v in q]
    print(f'kept faces: rest z {min(zs):.1f} to {max(zs):.1f}; x<0 {sum(1 for q, _ in keep if X[q[0]][0] < 0)}, '
          f'x>0 {sum(1 for q, _ in keep if X[q[0]][0] > 0)}; new nodes in the contact '
          f'{len(set(v for q, _ in keep for v in q) - set(v for q in allq[:len(have)] for v in q))}')
    if dry:
        return
    n0 = len(have)
    rows = ''.join(f'\t\t\t<quad4 id="{n0 + i}">{",".join(map(str, q))}</quad4>\n' for i, (q, _) in enumerate(keep, 1))
    t = t[:ms.end() - len('\t\t</Surface>\n')] + rows + t[ms.end() - len('\t\t</Surface>\n'):]
    msg = (f'NOT IN SOURCE: curve_all (curve_all_self) + {len(keep)} partly-tied curve boundary faces with every node at rest '
           f'z < {zcut:g} (of {len(cand)}; {len(drop)} left out: a node within {sr:g} mm of a non-neighbouring curve face at rest), '
           f'wound outward from their hex (tools/build_batch361.py); contact settings unchanged')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch361.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: float(a[a.index(k_) + 1]) if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--zcut', -12.0), g('--sr', 0.5), '--dry' in a)
