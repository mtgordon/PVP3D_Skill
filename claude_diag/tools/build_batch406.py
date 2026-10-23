"""Batch 406 (2026-10-22; the user: "I wanted the distal end to be filleted where it opens" + a FEBio Studio view looking into the
canal's distal opening: the edge where the curves' LUMEN-side surface (the canal opening) meets their distal end face).
A radius-R fillet of that edge, by node moves only (NOT IN SOURCE; connectivity, ties and surfaces unchanged):
  - lumen faces = the curves' boundary quads in the canal lumen contact (SlidingElastic1Primary / Secondary); end faces = the
    boundary quads within 1.5 mm of the end along the end tangent t with their normal along it; the opening rim = the end-face
    edges shared with a lumen face; n = the lumen faces' outward normal (into the lumen) at the rim;
  - a lumen-surface node at axial distance s < R from the end moves by -n (R - sqrt(R^2 - (R - s)^2));
  - an end-face node at depth d < R from the lumen surface (along -n) moves by -t (R - sqrt(R^2 - (R - d)^2));
  - a rim node (s ~ 0, d ~ 0) moves along -(t + n) / |t + n| by (sqrt(2) - 1) R (the arc's midpoint);
  so the corner becomes a quarter circle of radius R tangent to both surfaces. Only nodes of the curves alone move (nodes
  shared with a wall / the strap, or in a tie, keep their place); each move is cut back until every element containing the
  node keeps >= QMIN of its rest volume.
usage: py -3.10 build_batch406.py SRC_RUN NEW_RUN [--r 2.0] [--qmin 0.6] [--morph (interior curve nodes follow, harmonic)] [--holdouter (the end face's outer rim stays)] [--dry]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from build_batch297 import HF  # noqa: E402
from build_batch405 import vol, T5, W3  # noqa: E402
from junction_view_data import load  # noqa: E402

CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')


def unit(v):
    return v / np.linalg.norm(v)


def main(src, new, R, qmin, dry, morph=False, holdouter=False):
    t, X0, E = load(src)
    X = {k: v.copy() for k, v in X0.items()}
    hx = [c for nm, c in E.values() if nm in CURVES and len(c) == 8]
    cnt = Counter(frozenset(c[i] for i in f) for c in hx for f in HF)
    bq = [[c[i] for i in f] for c in hx for f in HF if cnt[frozenset(c[i] for i in f)] == 1]
    lumen = set()
    for nm in ('SlidingElastic1Primary', 'SlidingElastic1Secondary'):
        a = t.index(f'<Surface name="{nm}">'); b = t.index('</Surface>', a)
        lumen |= {frozenset(int(v) for v in s.split(',')) for s in re.findall(r'>([0-9,]+)</', t[a:b])}
    dn = defaultdict(set)
    for nm, c in E.values():
        for v in c:
            dn[v].add(nm)
    c0, c1 = t.find('<Constraints>'), t.find('</Constraints>')
    cons = {int(v) for v in re.findall(r'<node id="(\d+)"', t[c0:c1])}
    free = lambda v: not (dn[v] - set(CURVES)) and v not in cons   # noqa: E731
    owner = defaultdict(list)
    for k, (nm, c) in E.items():
        for v in c:
            owner[v].append(k)
    nrm = lambda q: unit(np.cross(X0[q[2]] - X0[q[0]], X0[q[3]] - X0[q[1]]))   # noqa: E731
    target, info = {}, []
    for sg in (-1, 1):
        cn = np.array([X0[v] for c in hx for v in c if np.sign(X0[v][0]) == sg])
        lo = cn[:, 2].min(); e = cn[cn[:, 2] < lo + 3].mean(0); u = cn[(cn[:, 2] > lo + 8) & (cn[:, 2] < lo + 12)].mean(0)
        tan = unit(e - u)
        smax = max(X0[v] @ tan for c in hx for v in c if np.sign(X0[v][0]) == sg)
        end = [q for q in bq if np.sign(X0[q[0]][0]) == sg and all(smax - X0[v] @ tan < 1.5 for v in q) and nrm(q) @ tan > 0.8]
        lum = [q for q in bq if frozenset(q) in lumen and np.sign(X0[q[0]][0]) == sg]
        endn = {v for q in end for v in q}
        lumn = {v for q in lum for v in q}
        rim = endn & lumn
        if not rim:
            continue
        tE = unit(np.mean([nrm(q) for q in end], 0))
        nL = unit(np.mean([nrm(q) for q in lum if set(q) & rim], 0))
        rimP = np.array([X0[v] for v in rim])
        for v in rim:
            if free(v):
                target[v] = X0[v] - (np.sqrt(2) - 1) * R * unit(tE + nL)
        for v in lumn - rim:
            s = smax - X0[v] @ tan
            # (L408 check) a lumen node on the end ring that is not a rim node sits beside the junction-tied end faces (their
            # end faces are not 'end' faces here): the formula would drop it ~R and fold its hexes -> it stays
            if 0.25 * R <= s < R and free(v):
                target[v] = X0[v] - nL * (R - np.sqrt(R * R - (R - s) ** 2))
        # (L408, the user: 'the outer edge doesn't need to be filleted in') the end face's OUTER rim stays: end-face moves taper
        # to 0 toward it (w = distance to the outer rim / R, capped at 1)
        ecnt = Counter(frozenset((q[k], q[(k + 1) % 4])) for q in end for k in range(4))
        outer = {v for e_, n_ in ecnt.items() if n_ == 1 for v in e_} - rim
        outP = np.array([X0[v] for v in outer]) if outer else None
        for v in endn - rim:
            d = min((p - X0[v]) @ nL for p in rimP)
            if 0 <= d < R and free(v):
                w = 1.0 if (not holdouter or outP is None) else min(1.0, np.min(np.linalg.norm(outP - X0[v], axis=1)) / R)
                target[v] = X0[v] - w * tE * (R - np.sqrt(R * R - (R - d) ** 2))
        info.append(f'{"left" if sg < 0 else "right"}: {len(end)} end faces, {len(lum)} lumen faces, rim {len(rim)} nodes '
                    f'({sum(free(v) for v in rim)} free), t {np.round(tE, 2)}, n {np.round(nL, 2)}')
    rest = {k: vol(X0, c) for k, (nm, c) in E.items()}
    if morph:
        # harmonic morph (L406, the fillet is deeper than the curve wall at the opening): surface nodes of the curves get their
        # fillet displacement (0 if not targeted), tied / shared nodes and nodes beyond 3R of the rim 0, interior curve nodes the
        # average of their hex-edge neighbours (Jacobi), so the elements behind the rounded edge follow it
        surf = {v for q in bq for v in q}
        rimall = np.array([X0[v] for v in target]) if target else np.zeros((1, 3))
        cnodes = {v for c in hx for v in c}
        zone = {v for v in cnodes if np.min(np.linalg.norm(rimall - X0[v], axis=1)) < 3 * R}
        nbr = defaultdict(set)
        for c in hx:
            for a_, b_ in ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)):
                nbr[c[a_]].add(c[b_]); nbr[c[b_]].add(c[a_])
        U = {v: (target[v] - X0[v]) if v in target else np.zeros(3) for v in cnodes}
        inner = [v for v in zone if v not in surf and free(v)]
        for _ in range(400):
            for v in inner:
                U[v] = np.mean([U[w] for w in nbr[v]], 0)
        target = {v: X0[v] + U[v] for v in zone if np.linalg.norm(U[v]) > 1e-6}
    moved, cut = {}, 0
    for v in sorted(target, key=lambda v: -np.linalg.norm(target[v] - X0[v])):
        ring2 = {k2 for k in owner[v] for w in E[k][1] for k2 in owner[w]} - set(owner[v])

        def ok(s_):
            X[v] = X0[v] + s_ * (target[v] - X0[v])
            if not all(rest[k] is None or vol(X, E[k][1]) >= qmin * rest[k] for k in owner[v]):
                return False
            for k in ring2:   # (L408 check) the node must not end inside a neighbouring element (a fold the volume misses)
                c = E[k][1]
                sp = T5 if len(c) == 8 else (W3 if len(c) == 6 else None)
                for q in sp or ():
                    a_, b_, c_, d_ = (X[c[i]] for i in q)
                    lam = np.linalg.solve(np.column_stack([b_ - a_, c_ - a_, d_ - a_]), X[v] - a_)
                    if min(1 - lam.sum(), lam.min()) > 1e-3:
                        return False
            return True
        if ok(1.0):
            s_ = 1.0
        else:
            lo_, hi_ = 0.0, 1.0
            for _ in range(20):
                m_ = (lo_ + hi_) / 2
                lo_, hi_ = (m_, hi_) if ok(m_) else (lo_, m_)
            s_ = lo_; cut += 1
        X[v] = X0[v] + s_ * (target[v] - X0[v])
        if np.linalg.norm(X[v] - X0[v]) > 1e-6:
            moved[v] = np.linalg.norm(X[v] - X0[v])
    # (L408 check) final pass: a node moved early can end inside a neighbour moved later; halve both moves until no moved
    # node lies inside an element it does not belong to
    from scipy.spatial import cKDTree
    for it in range(12):
        tets, tk = [], []
        for k, (nm, c) in E.items():
            sp = T5 if len(c) == 8 else (W3 if len(c) == 6 else None)
            if sp and set(c) & set(moved):
                tets += [[c[i] for i in s] for s in sp]; tk += [k] * len(sp)
        if not tets:
            break
        TP = np.array([[X[w] for w in q] for q in tets]); cen = TP.mean(1); tr = cKDTree(cen)
        RM = np.max(np.linalg.norm(TP - cen[:, None], axis=2))
        bad = set()
        for v in moved:
            for j in tr.query_ball_point(X[v], RM):
                if v in tets[j]:
                    continue
                a_, b_, c_, d_ = TP[j]
                lam = np.linalg.solve(np.column_stack([b_ - a_, c_ - a_, d_ - a_]), X[v] - a_)
                if min(1 - lam.sum(), lam.min()) > 1e-3:
                    bad.add(v); bad |= {w for w in E[tk[j]][1] if w in moved}
                    break
        # (L408 check) and no curve_all node within 0.5 mm BEHIND a curve_all face (the self-contact would grab it at rest)
        a0 = t.index('<Surface name="curve_all">'); b0 = t.index('</Surface>', a0)
        CA = [[int(w) for w in s.split(',')] for s in re.findall(r'>([0-9,]+)</', t[a0:b0])]
        can = sorted({w for f in CA for w in f}); ctr = cKDTree(np.array([X[w] for w in can]))
        for f in CA:
            if not set(f) & set(moved):
                continue
            Q = np.array([X[w] for w in f])
            for tri in ((0, 1, 2),) if len(f) == 3 else ((0, 1, 2), (0, 2, 3)):
                A, B, C = Q[list(tri)]; n = np.cross(B - A, C - A); n /= np.linalg.norm(n)
                for j in ctr.query_ball_point(Q.mean(0), 3.0):
                    w = can[j]
                    if w in f:
                        continue
                    d = (X[w] - A) @ n
                    if -0.5 < d < 0:
                        q = X[w] - d * n; v0, v1, v2 = B - A, C - A, q - A
                        d00, d01, d11, d20, d21 = v0 @ v0, v0 @ v1, v1 @ v1, v2 @ v0, v2 @ v1
                        den = d00 * d11 - d01 * d01
                        l1 = (d11 * d20 - d01 * d21) / den; l2 = (d00 * d21 - d01 * d20) / den
                        if min(1 - l1 - l2, l1, l2) > -0.02:
                            bad |= {w2 for w2 in list(f) + [w] if w2 in moved}
        if not bad:
            break
        for v in bad:
            X[v] = X0[v] + 0.5 * (X[v] - X0[v]); moved[v] = np.linalg.norm(X[v] - X0[v])
        cut += len(bad)
    mv = np.array(list(moved.values()) or [0])
    print('; '.join(info))
    print(f'opening fillet R {R:g}: {len(target)} nodes targeted, {len(moved)} moved (median {np.median(mv):.2f}, max {mv.max():.2f} mm), '
          f'{cut} cut back to keep every element >= {qmin:g} of its rest volume')
    if dry:
        return
    def node_sub(m):
        k = int(m.group(1))
        return f'{m.group(0)[:m.start(2) - m.start(0)]}{X[k][0]:.7g},{X[k][1]:.7g},{X[k][2]:.7g}</node>' if k in moved else m.group(0)
    t = re.sub(r'<node id="(\d+)">([^<]+)</node>', node_sub, t)
    msg = (f'NOT IN SOURCE (the user 10-22): the curves\' distal OPENING rounded (the edge between their lumen surface and their '
           f'end face) with a radius-{R:g} mm quarter-circle fillet: {len(moved)} free curve nodes moved (median {np.median(mv):.2f}, '
           f'max {mv.max():.2f} mm){" incl. interior nodes (harmonic morph within 3R)" if morph else ""}, {cut} cut back to keep every element >= {qmin:g} of its rest volume; ' + '; '.join(info) +
           '; node positions only')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch406.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--r', 2.0), g('--qmin', 0.6), '--dry' in a, '--morph' in a, '--holdouter' in a)
