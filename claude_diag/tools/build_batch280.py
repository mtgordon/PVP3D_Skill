"""Batch 280 (2026-10-11; the user: "get rid of the flares and make the curves attach better to the PeB ... try some things out").
On the rebuilt junction (L278_junction_*, build_batch278), MAX showed why the curve end still balloons: the wall everts past
the hymen and rolls over its attachment to the PeB (PVW / curves above move 16-19 mm, the tied curve end 13 mm). The thin PVW
bends through it; the 2.3 mm curve U takes the same roll-over with 2-3 coarse hex rings, which are crushed along the canal
(edges to 36 %) and bulge sideways (to 2.3x). Here (NOT IN SOURCE, a mesh refinement only): each of the curves' last NR PVW-side
rings (both sides) is split into NS thinner rings along the canal (new nodes on the axial edges, linear), so the bend is spread
over NR x NS rings. Split hexes keep their id for the first sub-ring; the others go into a new domain canal_wrap_hex_pvw_fine
(same material Curve_soft, same mass damping). Surface facets on split hexes are split too (every surface renumbered). New
nodes on edges whose ends are both attached to the PVW (PVW nodes, or curve nodes tied onto the PVW) are tied onto the PVW
face they lie on (junction_ties); new nodes on edges whose ends are both tied onto the PeB likewise onto the PeB.
usage: py -3.10 build_batch280.py SRC_RUN NEW_RUN [--rings NR | --zmax Z] [--split NS]
(--zmax Z: every whole curve ring, AVW and PVW halves, with a node shallower than Z mm, instead of the last NR PVW-side rings)"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load, frame  # noqa: E402
from junction_proposal import mirror_ids  # noqa: E402
from build_batch278 import bilinear, lc, RIGHT, PVW, PEB  # noqa: E402
from apex_close import HF  # noqa: E402
from build_batch162 import write  # noqa: E402

RIGHT_UP = [(27783, 27794)]          # the ring above ring A (right side), used when NR = 4


def rings_below(E, L, zmax):
    """every whole curve ring (canal_wrap_hex_avw + _pvw hexes joined through their lateral faces, axial edges i - i+4) with a
    node shallower than zmax"""
    cur = {k: c for k, (n, c) in E.items() if n in ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')}
    LAT = [(0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    byface = defaultdict(list)
    for k, c in cur.items():
        for f in LAT:
            byface[frozenset(c[i] for i in f)].append(k)
    seen, out = set(), []
    for k in cur:
        if k in seen:
            continue
        ring, stack = [], [k]
        seen.add(k)
        while stack:
            q = stack.pop(); ring.append(q)
            for f in LAT:
                for r in byface[frozenset(cur[q][i] for i in f)]:
                    if r not in seen:
                        seen.add(r); stack.append(r)
        if min(L[v][2] for q in ring for v in cur[q]) < zmax:
            out += ring
    return out


def main(src, new, nr=3, ns=3, zmax=None):
    t, X, E = load(src)
    c0, B = frame(X, E)
    L = {k: (p - c0) @ B for k, p in X.items()}
    if zmax is None:
        rings = (RIGHT_UP + RIGHT)[-nr:] if nr <= 4 else RIGHT
        ids = [k for a, b in rings for k in range(a, b + 1)]
        ids += mirror_ids(E, L, ids)
        assert all(E[k][0] == 'canal_wrap_hex_pvw' for k in ids)
    else:
        ids = rings_below(E, L, zmax)
        ax = [abs(((L[E[k][1][i + 4]] - L[E[k][1][i]]) / np.linalg.norm(L[E[k][1][i + 4]] - L[E[k][1][i]]))[2])
              for k in ids for i in range(4)]
        print(f'{len(ids)} hexes in whole rings below depth {zmax}; axial edges |cos to the canal axis| min {min(ax):.2f}, '
              f'median {np.median(ax):.2f}')
    dn = lambda d: {v for n, c in E.values() if n == d for v in c}   # noqa: E731
    pvw, peb, avw = dn(PVW), dn(PEB), dn('_PickedSet347')
    SETS = {'PVW': pvw, 'PeB': peb, 'AVW': avw}
    # nodes held on the PVW / on the PeB by junction_ties or canal_wrap_ends (by their host nodes)
    held = {}
    for name in ('junction_ties', 'canal_wrap_ends'):
        blk = re.search(r'<constraint name="%s".*?</constraint>' % name, t, re.S)
        for lcs in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', blk.group(0), re.S):
            nn = [int(a) for a in re.findall(r'<node id="(\d+)" bc="x">', lcs)]
            if nn:
                hosts = set(nn[1:])
                held[nn[0]] = next((p for p, S_ in SETS.items() if hosts <= S_), None)
    on = lambda v, part: v in SETS[part] or held.get(v) == part   # noqa: E731
    nid = max(X) + 1
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    newnode, newpos = {}, {}

    def mid(a, b, j):
        nonlocal nid
        key = (a, b, j) if a < b else (b, a, ns - j)
        if key not in newnode:
            f = key[2] / ns
            newnode[key] = nid
            newpos[nid] = (1 - f) * L[key[0]] + f * L[key[1]]
            nid += 1
        return newnode[key]
    layers = {}                          # element -> list of ns + 1 node layers (4 nodes each)
    for k in ids:
        c = E[k][1]
        lay = [c[:4]] + [[mid(c[i], c[i + 4], j) for i in range(4)] for j in range(1, ns)] + [c[4:]]
        layers[k] = lay
    sub = {}                             # element -> [(id, conn)]
    for k, lay in layers.items():
        sub[k] = [(k if j == 0 else None, lay[j] + lay[j + 1]) for j in range(ns)]
    fine = []
    for k in ids:
        for j, (i, c) in enumerate(sub[k]):
            if i is None:
                sub[k][j] = (eid, c); fine.append((eid, c)); eid += 1
    # ties for the new nodes
    ties, nt = '', defaultdict(int)
    pfaces = {}
    for part, dom in (('PVW', PVW), ('PeB', PEB), ('AVW', '_PickedSet347')):
        els = [c for n, c in E.values() if n == dom]
        cnt = defaultdict(int)
        for c in els:
            for fc in HF:
                cnt[frozenset(c[i] for i in fc)] += 1
        pfaces[part] = [[c[i] for i in fc] for c in els for fc in HF if cnt[frozenset(c[i] for i in fc)] == 1
                        or (part == 'PeB' and set(c[i] for i in fc) <= (pvw & peb))]
    worst = 0.0
    for (a, b, j), v in newnode.items():
        for part in ('PVW', 'PeB', 'AVW'):
            if on(a, part) and on(b, part):
                F = pfaces[part]
                C = np.array([np.mean([L[q] for q in f], 0) for f in F])
                best = None
                for i in np.argsort(np.linalg.norm(C - newpos[v], axis=1))[:30]:
                    Q = np.array([L[q] for q in F[i]])
                    w, d, m = bilinear(Q, newpos[v])
                    if m <= 1.0 + 1e-3 and (best is None or d < best[1]):
                        best = (list(zip(F[i], w)), d, w @ Q)
                if best is None:
                    continue
                worst = max(worst, best[1])
                newpos[v] = best[2]
                ties += lc(v, best[0]); nt[part] += 1
                break
    # ---- write: nodes, elements, domain, damping, surfaces, ties
    Xg = {v: c0 + B @ p for v, p in newpos.items()}
    a = t.index('<Nodes name="junction_pvw_bottom">'); b = t.index('</Nodes>', a) + len('</Nodes>\n')
    t = t[:b] + '\t\t<Nodes name="curve_end_fine_nodes">\n' + ''.join(
        f'\t\t\t<node id="{v}">' + ','.join(f'{c:.7g}' for c in Xg[v]) + '</node>\n' for v in sorted(Xg)) + '\t\t</Nodes>\n' + t[b:]
    for k in ids:
        old = f'<elem id="{k}">' + ','.join(map(str, E[k][1])) + '</elem>'
        assert old in t, k
        t = t.replace(old, f'<elem id="{k}">' + ','.join(map(str, sub[k][0][1])) + '</elem>', 1)
    a = t.index('<Elements type="hex8" name="apex_contour">'); b = t.index('</Elements>', a) + len('</Elements>\n')
    t = t[:b] + '\t\t<Elements type="hex8" name="canal_wrap_hex_pvw_fine">\n' + ''.join(
        f'\t\t\t<elem id="{i}">' + ','.join(map(str, c)) + '</elem>\n' for i, c in fine) + '\t\t</Elements>\n' + t[b:]
    t = t.replace('\t\t<SolidDomain name="canal_wrap_hex_pvw" mat="Curve_soft" />\n',
                  '\t\t<SolidDomain name="canal_wrap_hex_pvw" mat="Curve_soft" />\n'
                  '\t\t<SolidDomain name="canal_wrap_hex_pvw_fine" mat="Curve_soft" />\n', 1)
    m = re.search(r'\t\t<body_load name="wall_damping_canal_wrap_hex_pvw" .*?</body_load>\n', t, re.S)
    t = t[:m.end()] + m.group(0).replace('canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine') + t[m.end():]
    # surfaces: a facet with 2 nodes in a split hex's layer 0 and 2 in its layer ns -> ns facets
    lat = {}
    for k, lay in layers.items():
        pos = {}
        for i in range(4):
            pos[lay[0][i]] = (0, i); pos[lay[ns][i]] = (ns, i)
        for fc in HF:
            f = [(E[k][1])[q] for q in fc]
            if len({pos[v][0] for v in f}) == 2:
                lat[frozenset(f)] = (lay, pos)
    nsplit = 0

    def surf(m):
        nonlocal nsplit
        out, n = [], 0
        for fm in re.finditer(r'<(quad4|tri3) id="\d+">([\d,]+)</\1>', m.group(2)):
            f = [int(v) for v in fm.group(2).split(',')]
            key = frozenset(f)
            if fm.group(1) == 'quad4' and key in lat:
                lay, pos = lat[key]
                for j in range(ns):
                    g = [lay[j if pos[v][0] == 0 else j + 1][pos[v][1]] for v in f]
                    n += 1; out.append(f'\t\t\t<quad4 id="{n}">' + ','.join(map(str, g)) + '</quad4>\n')
                nsplit += 1
            else:
                n += 1; out.append(f'\t\t\t<{fm.group(1)} id="{n}">{fm.group(2)}</{fm.group(1)}>\n')
        return m.group(1) + ''.join(out) + '\t\t</Surface>'
    t = re.sub(r'(<Surface name="[^"]+">\n)(.*?)\t\t</Surface>', surf, t, flags=re.S)
    blk = re.search(r'\t\t<constraint name="junction_ties".*?</constraint>\n', t, re.S).group(0)
    t = t.replace(blk, blk.replace('\t\t</constraint>\n', ties + '\t\t</constraint>\n'))
    # check: smallest sub-tet of all the new hexes
    from curve_swell import TETS
    P = lambda v: newpos.get(v, L.get(v))   # noqa: E731
    mn = min(np.dot(P(c[q1]) - P(c[q0]), np.cross(P(c[q2]) - P(c[q0]), P(c[q3]) - P(c[q0]))) / 6
             for k in ids for _, c in sub[k] for q0, q1, q2, q3 in TETS)
    what = (f'the curves\' last {nr} PVW-side rings on each side' if zmax is None else
            f'every whole curve ring (AVW + PVW halves) reaching above depth {zmax} mm')
    msg = (f'NOT IN SOURCE (mesh refinement): {what} ({len(ids)} hexes) split '
           f'into {ns} rings each along the canal -> {len(ids) * ns} hexes ({len(fine)} new, domain canal_wrap_hex_pvw_fine, '
           f'Curve_soft); {len(newnode)} new nodes, {dict(nt)} of them tied onto the wall face they lie on (moved <= '
           f'{worst:.3f} mm); {nsplit} surface facets split; smallest sub-tet {mn:.3f} mm3')
    assert mn > 0
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch280.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    nr = int(a[a.index('--rings') + 1]) if '--rings' in a else 3
    ns = int(a[a.index('--split') + 1]) if '--split' in a else 3
    zmax = float(a[a.index('--zmax') + 1]) if '--zmax' in a else None
    main(a[0], a[1], nr, ns, zmax)
