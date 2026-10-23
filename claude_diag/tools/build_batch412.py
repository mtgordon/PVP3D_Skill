"""Batch 412 (2026-10-22; the user on L403: "for the corner edges between the avw and the pvw, there is a finer mesh towards the
distal end. Make that same mesh size for the entire corner edges." -> "build the whole-curve refinement on L403 so I can look
at it before it runs").
NOT IN SOURCE (mesh refinement only; the curves are ours): every whole curve ring (canal_wrap_hex_avw / _avw_cervix / _pvw /
_pvw_fine hexes joined through their lateral faces; axial edges i - i+4, as build_batch280) whose median axial edge is
> MAX mm (default 1.1) is split into ns = ceil(edge / MAX) thinner rings along the canal (the distal rings, ~1.0 mm from
build_batch280, stay). Through-wall layers and the around-the-U spacing are unchanged. The first sub-ring keeps the element
id; the others go into a new domain <parent>_len (the parent's material and mass damping).
New nodes: on an axial edge whose two ends are on / held on the PVW, PeB or AVW -> tied onto that part's face (as
build_batch280; junction_ties); else on an edge also used by an element that is not split (a wall element sharing curve nodes,
anything else) -> held at its point on the edge (linear constraints, curve_length_ties); else, when both ends are in a NodeSet,
added to it; otherwise free. Surface facets on split hex faces are split (every surface: contacts, the curve pressure, ...).
usage: py -3.10 build_batch412.py SRC_RUN NEW_RUN [--max 1.1] [--zmax Z] [--dry]
(--zmax Z: only rings with their centre at rest z <= Z; 10-22 L413 crawled with the whole curve split)"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch278 import bilinear, lc, PVW, PEB, EDGES  # noqa: E402
from apex_close import HF  # noqa: E402
from build_batch162 import write  # noqa: E402
from curve_swell import TETS  # noqa: E402

DOMS = ('canal_wrap_hex_avw', 'canal_wrap_hex_avw_cervix', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
LAT = [(0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def rings(E):
    cur = {k: c for k, (n, c) in E.items() if n in DOMS}
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
        out.append(ring)
    return out


def main(src, new, target=1.0, smin=1.1, dry=False, zmax=None):
    t, X, E = load(src)
    X = {k: np.asarray(p, float) for k, p in dict(X).items()}
    NS, ids = {}, []
    rep = []
    for ring in rings(E):
        ax = [np.linalg.norm(X[E[k][1][i + 4]] - X[E[k][1][i]]) for k in ring for i in range(4)]
        med = float(np.median(ax))
        zc = float(np.mean([X[v][2] for k in ring for v in E[k][1]]))
        xc = float(np.mean([X[v][0] for k in ring for v in E[k][1]]))
        ns = int(np.ceil(med / smin - 1e-9)) if med > smin else 1
        if zmax is not None and zc > zmax:
            ns = 1
        rep.append((xc > 0, zc, len(ring), med, ns, sorted({E[k][0] for k in ring})))
        if ns >= 2:
            for k in ring:
                NS[k] = ns
            ids += ring
    for side in (True, False):
        print('right (x > 0)' if side else 'left (x < 0)')
        for s, zc, n, med, ns, doms in sorted(r for r in rep if r[0] == side):
            print(f'  ring z {zc:6.1f}: {n:3d} hexes, axial median {med:.2f} mm -> {ns} ({",".join(d[11:] for d in doms)})')
    print(f'{len(ids)} hexes split, {sum(NS.values()) - len(ids)} new')
    if dry:
        return
    dn = lambda d: {v for n, c in E.values() if n == d for v in c}   # noqa: E731
    pvw, peb, avw = dn(PVW), dn(PEB), dn('_PickedSet347')
    SETS = {'PVW': pvw, 'PeB': peb, 'AVW': avw}
    held = {}
    for name in ('junction_ties', 'canal_wrap_ends'):
        blk = re.search(r'<constraint name="%s".*?</constraint>' % name, t, re.S)
        for lcs in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', blk.group(0), re.S):
            nn = [int(a) for a in re.findall(r'<node id="(\d+)" bc="x">', lcs)]
            if nn:
                held[nn[0]] = next((p for p, S_ in SETS.items() if set(nn[1:]) <= S_), None)
    on = lambda v, part: v in SETS[part] or held.get(v) == part   # noqa: E731
    split = set(ids)
    users = defaultdict(set)                      # edge -> elements using it (any element type)
    for k, (n, c) in E.items():
        ed = [(c[a], c[b]) for a, b in EDGES] if len(c) == 8 else [(c[i], c[(i + 1) % len(c)]) for i in range(len(c))]
        for a, b in ed:
            users[frozenset((a, b))].add(k)
    nsets = {m.group(1): {int(v) for v in re.findall(r'\d+', m.group(2))}
             for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', t)}
    nid = max(X) + 1
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    newnode, newpos, frac = {}, {}, {}

    def mid(a, b, j, ns):
        nonlocal nid
        key = (a, b, j, ns) if a < b else (b, a, ns - j, ns)
        if key not in newnode:
            f = key[2] / ns
            newnode[key] = nid; newpos[nid] = (1 - f) * X[key[0]] + f * X[key[1]]; frac[nid] = (key[0], key[1], f)
            nid += 1
        return newnode[key]
    layers = {}
    for k in ids:
        c, ns = E[k][1], NS[k]
        layers[k] = [c[:4]] + [[mid(c[i], c[i + 4], j, ns) for i in range(4)] for j in range(1, ns)] + [c[4:]]
    sub, fine = {}, defaultdict(list)
    for k, lay in layers.items():
        sub[k] = [(k, lay[0] + lay[1])]
        for j in range(1, NS[k]):
            sub[k].append((eid, lay[j] + lay[j + 1])); fine[E[k][0]].append((eid, lay[j] + lay[j + 1])); eid += 1
    # ties
    pfaces = {}
    for part, dom in (('PVW', PVW), ('PeB', PEB), ('AVW', '_PickedSet347')):
        els = [c for n, c in E.values() if n == dom]
        cnt = defaultdict(int)
        for c in els:
            for fc in HF:
                cnt[frozenset(c[i] for i in fc)] += 1
        pfaces[part] = [[c[i] for i in fc] for c in els for fc in HF if cnt[frozenset(c[i] for i in fc)] == 1
                        or (part == 'PeB' and set(c[i] for i in fc) <= (pvw & peb))]
    jties, eties, nt, worst, addset = '', '', defaultdict(int), 0.0, defaultdict(list)
    for v in sorted(frac):
        a, b, f = frac[v]
        done = False
        for part in ('PVW', 'PeB', 'AVW'):
            if on(a, part) and on(b, part):
                F = pfaces[part]
                C = np.array([np.mean([X[q] for q in fc], 0) for fc in F])
                best = None
                for i in np.argsort(np.linalg.norm(C - newpos[v], axis=1))[:30]:
                    Q = np.array([X[q] for q in F[i]])
                    w, d, m = bilinear(Q, newpos[v])
                    if m <= 1.0 + 1e-3 and (best is None or d < best[1]):
                        best = (list(zip(F[i], w)), d, w @ Q)
                if best is None:
                    continue
                worst = max(worst, best[1]); newpos[v] = best[2]
                jties += lc(v, best[0]); nt[part] += 1; done = True
                break
        if done:
            continue
        if users[frozenset((a, b))] - split:
            eties += lc(v, [(a, 1 - f), (b, f)]); nt['edge'] += 1
            continue
        hit = [s for s, S_ in nsets.items() if a in S_ and b in S_]
        for s in hit:
            addset[s].append(v)
        nt['nodeset' if hit else 'free'] += 1
    # check: smallest sub-tet of all the new hexes
    P = lambda v: newpos.get(v, X.get(v))   # noqa: E731
    mn = min(np.dot(P(c[q1]) - P(c[q0]), np.cross(P(c[q2]) - P(c[q0]), P(c[q3]) - P(c[q0]))) / 6
             for k in ids for _, c in sub[k] for q0, q1, q2, q3 in TETS)
    assert mn > 0, mn
    # ---- write: nodes after the last node block, elements, domains, damping, surfaces, node sets, ties
    b = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:b] + '\t\t<Nodes name="curve_length_nodes">\n' + ''.join(
        f'\t\t\t<node id="{v}">' + ','.join(f'{c:.7g}' for c in newpos[v]) + '</node>\n' for v in sorted(newpos)) + \
        '\t\t</Nodes>\n' + t[b:]
    for k in ids:
        old = f'<elem id="{k}">' + ','.join(map(str, E[k][1])) + '</elem>'
        assert t.count(old) == 1, k
        t = t.replace(old, f'<elem id="{k}">' + ','.join(map(str, sub[k][0][1])) + '</elem>', 1)
    me = t.index('</Mesh>'); b = t.rindex('</Elements>', 0, me) + len('</Elements>\n')
    t = t[:b] + ''.join(f'\t\t<Elements type="hex8" name="{d}_len">\n' + ''.join(
        f'\t\t\t<elem id="{i}">' + ','.join(map(str, c)) + '</elem>\n' for i, c in fine[d]) + '\t\t</Elements>\n'
        for d in DOMS if fine[d]) + t[b:]
    for d in DOMS:
        if not fine[d]:
            continue
        m = re.search(r'\t\t<SolidDomain name="%s" (mat="[^"]+") />\n' % d, t)
        t = t[:m.end()] + f'\t\t<SolidDomain name="{d}_len" {m.group(1)} />\n' + t[m.end():]
        m = re.search(r'\t\t<body_load name="wall_damping_%s" .*?</body_load>\n' % d, t, re.S)
        t = t[:m.end()] + m.group(0).replace(f'"wall_damping_{d}"', f'"wall_damping_{d}_len"').replace(
            f'elem_set="{d}"', f'elem_set="{d}_len"') + t[m.end():]
    lat = {}
    for k, lay in layers.items():
        ns = NS[k]
        pos = {}
        for i in range(4):
            pos[lay[0][i]] = (0, i); pos[lay[ns][i]] = (ns, i)
        for fc in HF:
            f = [E[k][1][q] for q in fc]
            if len({pos[v][0] for v in f}) == 2:
                lat[frozenset(f)] = (lay, pos, ns)
    nsplit = 0

    def surf(m):
        nonlocal nsplit
        out, n = [], 0
        for fm in re.finditer(r'<(quad4|tri3) id="\d+">([\d,]+)</\1>', m.group(2)):
            f = [int(v) for v in fm.group(2).split(',')]
            key = frozenset(f)
            if fm.group(1) == 'quad4' and key in lat:
                lay, pos, ns = lat[key]
                for j in range(ns):
                    g = [lay[j if pos[v][0] == 0 else j + 1][pos[v][1]] for v in f]
                    n += 1; out.append(f'\t\t\t<quad4 id="{n}">' + ','.join(map(str, g)) + '</quad4>\n')
                nsplit += 1
            else:
                n += 1; out.append(f'\t\t\t<{fm.group(1)} id="{n}">{fm.group(2)}</{fm.group(1)}>\n')
        return m.group(1) + ''.join(out) + '\t\t</Surface>'
    t = re.sub(r'(<Surface name="[^"]+">\n)(.*?)\t\t</Surface>', surf, t, flags=re.S)
    for s, vs in addset.items():
        m = re.search(r'(<NodeSet name="%s">)([^<]*)(</NodeSet>)' % re.escape(s), t)
        t = t[:m.start(2)] + m.group(2).rstrip() + ',' + ','.join(map(str, vs)) + t[m.end(2):]
    blk = re.search(r'\t\t<constraint name="junction_ties".*?</constraint>\n', t, re.S).group(0)
    t = t.replace(blk, blk.replace('\t\t</constraint>\n', jties + '\t\t</constraint>\n'))
    if eties:
        hdr = re.search(r'\t\t<constraint name="junction_ties"[^>]*>\n(.*?)(?=\t\t\t<linear_constraint>)', t, re.S)
        new_blk = (hdr.group(0).replace('name="junction_ties"', 'name="curve_length_ties"') + eties + '\t\t</constraint>\n')
        t = t.replace(blk.replace('\t\t</constraint>\n', jties + '\t\t</constraint>\n'),
                      blk.replace('\t\t</constraint>\n', jties + '\t\t</constraint>\n') + new_blk, 1)
    msg = (f'NOT IN SOURCE (mesh refinement): every whole curve ring with a median axial edge >= {smin} mm split along the '
           f'canal to ~{target} mm ({len(ids)} hexes -> {sum(NS[k] for k in ids)}; {sum(len(v) for v in fine.values())} new in '
           f'{", ".join(d + "_len" for d in DOMS if fine[d])}); {len(newpos)} new nodes: {dict(nt)} (wall-face ties moved '
           f'<= {worst:.3f} mm; NodeSets extended: {dict((s, len(v)) for s, v in addset.items())}); {nsplit} surface facets '
           f'split; smallest sub-tet {mn:.4f} mm3')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch412.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    tg = float(a[a.index('--target') + 1]) if '--target' in a else 1.0
    sm = float(a[a.index('--max') + 1]) if '--max' in a else (float(a[a.index('--min') + 1]) if '--min' in a else 1.1)
    zm = float(a[a.index('--zmax') + 1]) if '--zmax' in a else None
    main(a[0], a[1], tg, sm, '--dry' in a, zm)
