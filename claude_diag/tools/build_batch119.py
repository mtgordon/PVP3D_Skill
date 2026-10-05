"""Batch 119 (2026-09-29, the user's request): the canal as a tube, version 2: the rounded wrap (NOT IN SOURCE).
L118_tube_fold (the AVW / cervix and PVW joined along both sides of the canal with shared lumen-edge nodes) + a half-disk
of wall tissue on each side, so the walls wrap round like the cervix into the PVW at the top:
  * at every AVW / cervix seam station (the 58 lumen-edge nodes a side), the cross-section's half-disk is centred on the
    closed lumen corner c (the seam node); it spans from the AVW's side face (its node column c, a1, a2 through the wall)
    round laterally to the PVW's side face (its column c, p1, p2, p3); radius ~ the wall thickness;
  * rings: the union of the two walls' layer heights (as fractions of each wall's thickness: the AVW's a1 and the PVW's
    p1, p2, sorted) plus the outer surface, so every real wall node on the side faces is a wrap node: at the AVW end the
    rings at the AVW's own heights share its nodes, the others lie on its column and follow it (linear constraints between
    the two column nodes around them); the same at the PVW end; 6 segments round the half-disk; wedges (penta6) at the
    centre, hexes (hex8) outside;
  * along the canal the wrap uses the AVW's stations; at those between two merged PVW columns (the PVW side is half as
    fine) the wrap's PVW-end nodes follow the PVW's side face (bilinear constraints on the two neighbouring columns);
    stations beyond the first / last merged PVW column at a canal end are left unwrapped;
  * material: the walls' (Vagina_AVW, the faithful Yeoh); constraints penalty 10, maxaug 10 (as LA_truss_ties).
  L119_tube_wrap   L118_tube_fold + the rounded wrap
usage: py -3.10 build_batch119.py [--check]   (--check: build in memory, print the mesh checks, write nothing)"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch90 import seam_nodes, _bedges  # noqa: E402
from build_batch118 import fold, BASE  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
NSEG = 6
LABEL = 'NOT IN SOURCE (the canal as a tube: rounded wrap)'


def hex_jac_ok(P):
    """Signed volume check of a hex8 (node order 0-3 bottom, 4-7 top): the Jacobian at the centre > 0."""
    d1 = (P[1] + P[2] + P[6] + P[5]) - (P[0] + P[3] + P[7] + P[4])
    d2 = (P[3] + P[2] + P[6] + P[7]) - (P[0] + P[1] + P[5] + P[4])
    d3 = (P[4] + P[5] + P[6] + P[7]) - (P[0] + P[1] + P[2] + P[3])
    return np.dot(d1, np.cross(d2, d3))


def wedge_jac(P):
    """penta6 (0-2 bottom triangle, 3-5 top): the triangle normal (bottom) against the extrusion direction."""
    n = np.cross(P[1] - P[0], P[2] - P[0])
    e = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3
    return np.dot(n, e)


def wrap(model):
    path = model.src                      # the base before the fold: columns from its geometry
    f0 = Feb(path)
    seam, fP, X = seam_nodes(Model7(path))
    eP = sorted({v for e in _bedges(fP) for v in e})
    PE = np.array([X[v] for v in eP])
    cx, peb = set(f0.domain_nodes(tube.CX)), set(f0.domain_nodes(tube.PEB))
    pedge = sorted({eP[int(np.argmin(np.linalg.norm(PE - X[v], axis=1)))] for v in seam})
    pedge = [p for p in pedge if p not in cx and p not in peb]
    colA = tube.columns(f0, (tube.AVW, tube.CX), 'SlidingElastic1Secondary', seam)
    colP = tube.columns(f0, (tube.PVW, tube.PEB), 'SlidingElastic1Primary', pedge)
    rep, ties = fold(model)               # the fold first (its merges and in-between constraints)
    inv = {a: p for p, a in rep.items()}  # merged seam node -> the PVW edge node it replaced
    thA = np.median([np.linalg.norm(X[c[-1]] - X[c[0]]) for c in colA.values() if len(c) == 3])
    thP = np.median([np.linalg.norm(X[c[-1]] - X[c[0]]) for c in colP.values() if len(c) == 4])

    def trim(col, th, n):
        """the column's first n nodes (the wall's layers: AVW 3, PVW 4), if it has them and they span <= 2.5 x the
        wall thickness (a walk that ran on along the rim or into the perineal body is cut to its first n nodes)"""
        if len(col) < n:
            return None
        c = col[:n]
        return c if np.linalg.norm(X[c[-1]] - X[c[0]]) <= 2.5 * th else None
    new_nodes, cons, hexes, wedges = [], [], [], []
    nid0 = model.max_node_id() + 1

    def add(p):
        new_nodes.append(p)
        return nid0 + len(new_nodes) - 1
    pos = dict(X)
    report = {}
    for sgn in (-1, 1):
        S = [a for a in seam if np.sign(X[a][0]) == sgn]
        ax = np.linalg.svd(np.array([X[a] for a in S]) - np.mean([X[a] for a in S], axis=0))[2][0]
        S.sort(key=lambda a: X[a] @ ax)
        merged = [i for i, a in enumerate(S) if a in inv]
        if len(merged) < 2:
            continue
        stations = []
        for i in range(merged[0], merged[-1] + 1):
            a0 = S[i]
            A = trim(colA.get(a0, []), thA, 3)
            if A is None:
                stations.append(None)
                continue
            # the PVW column: this station's own (merged) or the two merged neighbours' (for bilinear constraints)
            if a0 in inv and inv[a0] in colP:
                Pc = trim(colP[inv[a0]], thP, 4)
                pv = ('own', Pc)
            else:
                lo = max(j for j in merged if j < i)
                hi = min(j for j in merged if j > i)
                Pl, Ph = trim(colP.get(inv[S[lo]], []), thP, 4), trim(colP.get(inv[S[hi]], []), thP, 4)
                if Pl is None or Ph is None:
                    stations.append(None)
                    continue
                s = (X[a0] - X[S[lo]]) @ ax / max(1e-9, (X[S[hi]] - X[S[lo]]) @ ax)
                s = float(np.clip(s, 0, 1))
                pv = ('between', Pl, Ph, s)
            if pv[0] == 'own' and pv[1] is None:
                stations.append(None)
                continue
            stations.append((a0, A, pv))
        # build the rings station by station
        rings = []                         # per station: node ids [ring][segment]
        for st in stations:
            if st is None:
                rings.append(None)
                continue
            a0, A, pv = st
            c = X[a0]
            if pv[0] == 'own':
                Pn = pv[1]
                Pp = [X[v] for v in Pn]
            else:
                _, Pl, Ph, s = pv
                Pp = [(1 - s) * X[u] + s * X[w] for u, w in zip(Pl, Ph)]
                Pp[0] = c
            up = X[A[-1]] - c
            dn = Pp[-1] - c
            RA, RP = np.linalg.norm(up), np.linalg.norm(dn)
            up, dn = up / RA, dn / RP
            fa = [np.linalg.norm(X[A[1]] - c) / RA]
            fp = [np.linalg.norm(Pp[1] - c) / RP, np.linalg.norm(Pp[2] - c) / RP]
            fr = sorted([(x, 'A1') for x in fa] + [(x, 'P1') for x in fp[:1]] + [(x, 'P2') for x in fp[1:]]) + [(1.0, 'O')]
            # lateral outward direction: perpendicular to the column, pointing away from the canal's midline
            t = ax - (ax @ up) * up
            w = np.cross(up, t)
            w = w - (w @ up) * up
            w /= np.linalg.norm(w)
            if w[0] * sgn < 0:
                w = -w
            ring_ids = []
            for k, (fk, kind) in enumerate(fr):
                seg_ids = []
                for j in range(NSEG + 1):
                    th = np.pi / 2 - np.pi * j / NSEG           # +90 (AVW end) .. -90 (PVW end)
                    vert = up if th >= 0 else dn
                    d = np.cos(th) * w + abs(np.sin(th)) * vert
                    d /= np.linalg.norm(d)
                    tj = j / NSEG
                    r = (1 - tj) * fk * RA + tj * fk * RP
                    if j == 0:                                   # the AVW end: share its nodes where they exist
                        if kind == 'A1':
                            seg_ids.append(A[1]); continue
                        if kind == 'O':
                            seg_ids.append(A[-1]); continue
                        p = c + fk * RA * up
                        v = add(p)
                        lo_n, hi_n = (A[0], A[1]) if fk < fa[0] else (A[1], A[2])
                        f0_, f1_ = (0.0, fa[0]) if fk < fa[0] else (fa[0], 1.0)
                        ss = (fk - f0_) / (f1_ - f0_)
                        cons.append((v, [(lo_n, 1 - ss), (hi_n, ss)]))
                        seg_ids.append(v); continue
                    if j == NSEG:                                # the PVW end
                        if pv[0] == 'own':
                            if kind in ('P1', 'P2', 'O'):
                                seg_ids.append({'P1': Pn[1], 'P2': Pn[2], 'O': Pn[3]}[kind]); continue
                            p = c + fk * RP * dn
                            v = add(p)
                            lo_n, hi_n = (Pn[1], Pn[2]) if fp[0] <= fk <= fp[1] else ((a0, Pn[1]) if fk < fp[0] else (Pn[2], Pn[3]))
                            ya, yb = (np.linalg.norm(X[lo_n] - c) if lo_n != a0 else 0.0), np.linalg.norm(X[hi_n] - c)
                            ss = float(np.clip((fk * RP - ya) / max(1e-9, yb - ya), 0, 1))
                            cons.append((v, [(lo_n, 1 - ss), (hi_n, ss)]))
                            seg_ids.append(v); continue
                        _, Pl, Ph, s = pv
                        p = c + fk * RP * dn
                        v = add(p)
                        # bilinear on the two neighbouring PVW columns at the same fraction of the thickness
                        wts = []
                        for col, wc in ((Pl, 1 - s), (Ph, s)):
                            yy = [np.linalg.norm(X[q] - X[col[0]]) for q in col]
                            y = fk * yy[-1]
                            m_ = int(np.clip(np.searchsorted(yy, y) - 1, 0, len(col) - 2))
                            u_ = float(np.clip((y - yy[m_]) / max(1e-9, yy[m_ + 1] - yy[m_]), 0, 1))
                            q0 = col[m_] if m_ > 0 else rep.get(col[0], col[0])
                            wts += [(q0, wc * (1 - u_)), (col[m_ + 1], wc * u_)]
                        cons.append((v, wts))
                        seg_ids.append(v); continue
                    seg_ids.append(add(c + r * d))
                ring_ids.append(seg_ids)
            rings.append((a0, ring_ids))
        # elements between consecutive wrapped stations
        nh = nw = 0
        allpos = lambda ids: np.array([pos[i] if i in pos else new_nodes[i - nid0] for i in ids])
        for i in range(len(rings) - 1):
            if rings[i] is None or rings[i + 1] is None:
                continue
            (c0_, R0), (c1_, R1) = rings[i], rings[i + 1]
            for j in range(NSEG):
                wd = [c0_, R0[0][j], R0[0][j + 1], c1_, R1[0][j], R1[0][j + 1]]
                if wedge_jac(allpos(wd)) < 0:
                    wd = [c0_, R0[0][j + 1], R0[0][j], c1_, R1[0][j + 1], R1[0][j]]
                wedges.append(wd); nw += 1
                for k in range(len(R0) - 1):
                    hx = [R0[k][j], R0[k][j + 1], R0[k + 1][j + 1], R0[k + 1][j],
                          R1[k][j], R1[k][j + 1], R1[k + 1][j + 1], R1[k + 1][j]]
                    if hex_jac_ok(allpos(hx)) < 0:
                        hx = hx[4:] + hx[:4]
                    hexes.append(hx); nh += 1
        report[sgn] = (sum(r is not None for r in rings), len(rings), nh, nw)
    # write: nodes, elements, domains, constraints
    ids = model.add_nodes('canal_wrap_nodes', new_nodes)
    assert ids[0] == nid0
    blocks = model.mesh.findall('Elements')
    start = max(int(e.get('id')) for b in blocks for e in b) + 1
    last = blocks[-1]
    eh = ET.Element('Elements', type='hex8', name='canal_wrap_hex')
    for q, h in enumerate(hexes):
        ET.SubElement(eh, 'elem', id=str(start + q)).text = ','.join(map(str, h))
    ew = ET.Element('Elements', type='penta6', name='canal_wrap_wedge')
    for q, wd in enumerate(wedges):
        ET.SubElement(ew, 'elem', id=str(start + len(hexes) + q)).text = ','.join(map(str, wd))
    i0 = list(model.mesh).index(last)
    model.mesh.insert(i0 + 1, eh)
    model.mesh.insert(i0 + 2, ew)
    doms = model.root.find('MeshDomains')
    for nm in ('canal_wrap_hex', 'canal_wrap_wedge'):
        ET.SubElement(doms, 'SolidDomain', name=nm, mat='Vagina_AVW')
    c = ET.SubElement(model.root.find('Constraints'), 'constraint', name='canal_wrap_ends', type='linear constraint')
    ET.SubElement(c, 'tol').text = '0.01'
    ET.SubElement(c, 'penalty').text = '10'
    ET.SubElement(c, 'maxaug').text = '10'
    for v, wts in cons:
        for dof in ('x', 'y', 'z'):
            lc = ET.SubElement(c, 'linear_constraint')
            ET.SubElement(lc, 'node', id=str(v), bc=dof).text = '1'
            for q, wq in wts:
                if abs(wq) > 1e-9:
                    ET.SubElement(lc, 'node', id=str(q), bc=dof).text = f'{-wq:.6f}'
    # checks
    allP = {**{k: v for k, v in X.items()}, **{nid0 + i: p for i, p in enumerate(new_nodes)}}
    jh = np.array([hex_jac_ok(np.array([allP[n] for n in h])) for h in hexes])
    jw = np.array([wedge_jac(np.array([allP[n] for n in w_])) for w_ in wedges])
    model.log.append(f'{LABEL}: on the fold: per side (stations wrapped / spanned, hexes, wedges) {report}; {len(new_nodes)} new '
                     f'nodes, {len(hexes)} hex8 + {len(wedges)} penta6 (domains canal_wrap_hex / canal_wrap_wedge, material '
                     f'Vagina_AVW); {len(cons)} wrap-end nodes following the wall side faces (constraint canal_wrap_ends); '
                     f'wall thickness at the sides AVW {thA:.2f} / PVW {thP:.2f} mm; {NSEG} segments round the half-disk; '
                     f'element centre Jacobians all positive: hex {bool((jh > 0).all())}, wedge {bool((jw > 0).all())}')
    return report, jh, jw


if __name__ == '__main__':
    name = 'L119_tube_wrap'
    m = Model7(os.path.join(RUNS, BASE, BASE + '.feb'))
    rep_, jh, jw = wrap(m)
    print(m.log[-1])
    if '--check' not in sys.argv:
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        emit(name, m)
