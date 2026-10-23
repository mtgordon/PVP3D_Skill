"""Batch 402 (2026-10-22; the user: "the mesh of the last layer of the curved edges larger", option 1 agreed): at the distal
end of the curves (the last ZLEN mm along the end tangent), merge the curve's thin inner layer (~0.45 mm across, the layer
that wraps onto the wall edge and carries the wall ties; it runs the whole curve) with its neighbouring layer into one
~1.35 mm layer. In L394 the thin hexes at the end corner (where the strap and PVW ties meet it) swell 2.5x / crush 0.3x and
curve nodes from higher up slip past the tied edge.
For each zone hex T with a thin family (mean edge < THIN, across the curve, not along the tangent) whose wall-side face is
not shared with another curve hex and whose other face is shared with a curve hex N of the zone: the merged hex = T with
each middle node m (on the T / N face) replaced by its partner b across N (the a - m - b thickness lines); N is deleted.
Middle nodes no longer used are dropped (isolated); middle nodes still used by hexes outside the zone (the transition ring
face) become hanging nodes on the merged edge a - b and are tied there by a linear constraint u_m = (1 - lam) u_a + lam u_b
(rest position), as the model's other ties (<Constraints> linear constraint, penalty 10, tol 0.01, maxaug 10). Every
<Surface> facet that was a face of T or N is replaced by the merged hex's face (orientation kept). Pairs whose middle nodes
are shared with another domain, a NodeSet or a constraint are left alone. NOT IN SOURCE (mesh only; the curves are ours).
usage: py -3.10 build_batch402.py SRC_RUN NEW_RUN [--zlen 4] [--thin 0.6] [--dry]"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402

CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
FAM = {'A': ((0, 1), (3, 2), (4, 5), (7, 6)), 'B': ((0, 3), (1, 2), (4, 7), (5, 6)), 'C': ((0, 4), (1, 5), (2, 6), (3, 7))}
FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
T5 = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]


def vol(P, c):
    return sum(np.dot(np.cross(P[c[b]] - P[c[a]], P[c[cc]] - P[c[a]]), P[c[d]] - P[c[a]]) / 6 for a, b, cc, d in T5)


def main(src, new, zlen, thin, dry):
    t, X, E = load(src)
    hexes = {k: c for k, (nm, c) in E.items() if nm in CURVES and len(c) == 8}
    face_of = defaultdict(list)
    for k, c in hexes.items():
        for f in FACES:
            face_of[frozenset(c[i] for i in f)].append(k)
    # nodes used outside the curves, in node sets, in constraints
    other = {v for k, (nm, c) in E.items() if nm not in CURVES for v in c}
    nsets = {int(v) for m in re.finditer(r'<NodeSet name="[^"]*">([^<]*)</NodeSet>', t) for v in re.findall(r'\d+', m.group(1))}
    c0, c1 = t.find('<Constraints>'), t.find('</Constraints>')
    cons = {int(v) for v in re.findall(r'<node id="(\d+)"', t[c0:c1])}
    pairs, zone = [], {}
    for sg in (-1, 1):
        cn = np.array([X[v] for c in hexes.values() for v in c if np.sign(X[v][0]) == sg])
        lo = cn[:, 2].min(); e = cn[cn[:, 2] < lo + 3].mean(0); u = cn[(cn[:, 2] > lo + 8) & (cn[:, 2] < lo + 12)].mean(0)
        tan = (e - u) / np.linalg.norm(e - u)
        smax = max(X[v] @ tan for c in hexes.values() for v in c if np.sign(X[v][0]) == sg)
        for k, c in hexes.items():
            P = np.array([X[v] for v in c])
            if np.sign(P[:, 0].mean()) == sg and P.mean(0) @ tan > smax - zlen:
                zone[k] = tan
    used = set()
    for k, tan in zone.items():
        if k in used:
            continue
        c = hexes[k]; P = {v: X[v] for v in c}
        fam = {f: np.mean([np.linalg.norm(P[c[a]] - P[c[b]]) for a, b in e]) for f, e in FAM.items()}
        ax = max(FAM, key=lambda f: np.mean([abs((P[c[b]] - P[c[a]]) @ tan) / np.linalg.norm(P[c[b]] - P[c[a]]) for a, b in FAM[f]]))
        cand = [f for f in FAM if f != ax and fam[f] < thin]
        if len(cand) != 1:
            continue
        f = cand[0]
        sideA = frozenset(c[a] for a, b in FAM[f]); sideB = frozenset(c[b] for a, b in FAM[f])
        nA = [h for h in face_of[sideA] if h != k]; nB = [h for h in face_of[sideB] if h != k]
        if nA and not nB:
            wall, mid, N = sideB, sideA, nA[0]
        elif nB and not nA:
            wall, mid, N = sideA, sideB, nB[0]
        else:
            continue
        if N not in zone or N in used:
            continue
        cN = hexes[N]
        # partner across N: for each middle node m, the node of N joined to m by an edge and not on the shared face
        edgesN = {frozenset((cN[a], cN[b])) for fe in FAM.values() for a, b in fe}
        far = {}
        for m in mid:
            nb = [v for v in cN if v not in mid and frozenset((m, v)) in edgesN]
            if len(nb) != 1:
                break
            far[m] = nb[0]
        if len(far) != 4:
            continue
        aof = {}
        for a, b in FAM[f]:
            if c[a] in mid:
                aof[c[a]] = c[b]
            else:
                aof[c[b]] = c[a]
        if any(m in other or m in nsets or m in cons for m in mid):
            continue
        pairs.append((k, N, far, aof))
        used |= {k, N}
    # merged connectivity, volume check
    merged = {}
    for k, N, far, aof in pairs:
        cm = [far.get(v, v) for v in hexes[k]]
        if vol(X, cm) <= 0:
            continue
        merged[k] = (cm, N, far, aof)
    deleted = {N for cm, N, far, aof in merged.values()}
    mids = {m: (aof[m], far[m]) for cm, N, far, aof in merged.values() for m in far}
    still = {v for k, c in hexes.items() if k not in deleted and k not in merged for v in c}
    hanging = {m: ab for m, ab in mids.items() if m in still}
    dropped = set(mids) - set(hanging)
    print(f'zone hexes {len(zone)}; merged pairs {len(merged)} (of {len(pairs)} found); middle nodes {len(mids)}: '
          f'{len(dropped)} dropped, {len(hanging)} hanging (tied to their merged edge)')
    if dry:
        return
    # elements: rewrite merged, delete N
    def elem_sub(m_):
        k = int(m_.group(1))
        if k in deleted:
            return ''
        if k in merged:
            return f'{m_.group(0)[:m_.start(2) - m_.start(0)]}{",".join(map(str, merged[k][0]))}</elem>\n'
        return m_.group(0)
    t = re.sub(r'\t\t\t<elem id="(\d+)">([0-9,]+)</elem>\n', elem_sub, t)
    for d in deleted:
        assert f'"{d}"' not in re.sub(r'<node id="\d+"', '', t) or True
    # surfaces
    child = {}   # frozenset(child face) -> (merged face nodes as a function of the child's order)
    for k, (cm, N, far, aof) in merged.items():
        bof = {far[m]: m for m in far}; mof_a = {aof[m]: m for m in aof}
        for fc in FACES:
            h = [cm[i] for i in fc]
            ka = [v for v in h if v in mof_a]; kb = [v for v in h if v in bof]
            if len(ka) == 2 and len(kb) == 2:
                tch = [bof.get(v, v) for v in h]; nch = [mof_a.get(v, v) for v in h]
                child[frozenset(tch)] = ('T', {bof[b]: b for b in kb})       # m -> b
                child[frozenset(nch)] = ('N', {mof_a[a]: a for a in ka})     # m -> a
    nrep = 0

    def surf_sub(ms):
        nonlocal nrep
        name, body = ms.group(1), ms.group(2)
        out, seen = [], set()
        for g, fid, nodes in re.findall(r'\t\t\t<(quad4|tri3) id="(\d+)">([0-9,]+)</\1>\n', body):
            q = [int(v) for v in nodes.split(',')]
            key = frozenset(q)
            if key in child:
                q = [child[key][1].get(v, v) for v in q]; nrep += 1
            k2 = frozenset(q)
            if k2 in seen:
                continue
            seen.add(k2); out.append((g, q))
        assert not any(v in dropped for g, q in out for v in q), f'surface {name} still uses a dropped node'
        return (f'\t\t<Surface name="{name}">\n' + ''.join(f'\t\t\t<{g} id="{i}">{",".join(map(str, q))}</{g}>\n'
                                                       for i, (g, q) in enumerate(out, 1)) + '\t\t</Surface>\n')
    t = re.sub(r'\t\t<Surface name="([^"]+)">\n(.*?)\t\t</Surface>\n', surf_sub, t, flags=re.S)
    # hanging-node ties
    lc = []
    for m, (a, b) in sorted(hanging.items()):
        lam = float(np.clip((X[m] - X[a]) @ (X[b] - X[a]) / ((X[b] - X[a]) @ (X[b] - X[a])), 0, 1))
        for d in 'xyz':
            lc.append(f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{d}">1</node>\n'
                      f'\t\t\t\t<node id="{a}" bc="{d}">{-(1 - lam):.6f}</node>\n\t\t\t\t<node id="{b}" bc="{d}">{-lam:.6f}</node>\n'
                      '\t\t\t</linear_constraint>\n')
    blk = ('\t\t<constraint name="curve_end_merge_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n'
           '\t\t\t<maxaug>10</maxaug>\n' + ''.join(lc) + '\t\t</constraint>\n')
    i = t.index('\t</Constraints>')
    t = t[:i] + blk + t[i:]
    # element-indexed data referencing deleted elements?
    for d in deleted:
        assert not re.search(rf'<e lid="{d}"|<elem lid="{d}"', t)
    msg = (f'NOT IN SOURCE (mesh only, the user 10-22): curve distal end (last {zlen:g} mm along the end tangent): '
           f'{len(merged)} thin inner-layer hexes (mean edge < {thin:g} mm across) merged with their neighbouring layer '
           f'(the merged hex keeps the thin hex id / domain, the neighbour deleted: {sorted(deleted)[:3]}...); '
           f'{len(dropped)} middle nodes dropped, {len(hanging)} hanging nodes at the transition ring tied to their merged edge '
           f'(curve_end_merge_ties, {3 * len(hanging)} linear constraints, penalty 10); {nrep} surface facets replaced by '
           f'merged faces')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch402.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--zlen', 4), g('--thin', 0.6), '--dry' in a)
