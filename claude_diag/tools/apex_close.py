"""Close the tube's apical corners (2026-10-07, the user: "it should be a completely closed tube at the apical side ... the
holes are on the corner where the current curve doesn't go all of the way to the back").
The curves (canal_wrap_hex_avw / _pvw) end at the last seam station (~-86 mm; 0 = hymen at rest); the cervix and PVW cap the
canal at ~-90..-93, so each lateral apical corner is open. Each curve is EXTRUDED beyond its apical end cap (24 faces,
35 nodes = 5 through-thickness x 7 round the U) along each cap node's own ring-to-ring direction, in K rings whose last one
lands on the cervix / PVW surface (each node's ray to the nearest cervix / PVW boundary face; fallback: its last ring step).
The new hexes are translated copies of the cap hexes (same node order, so the same orientation), in new domains
<curve>_apex with the curve half's material and damping. Every node of the last ring, and the end columns (nodes tied to /
shared with a wall at the cap) of every new ring, follow the nearest cervix / PVW node (linear constraints, as
canal_wrap_ends: penalty 10, maxaug 10). NOT IN SOURCE.
usage (module): t2, msg = close(t, K=3)"""
import re
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
CURVES = {'canal_wrap_hex_avw': 'Vagina_AVW', 'canal_wrap_hex_pvw': 'Vagina_PVW'}
CAPWALLS = ('_PickedSet346', '_PickedSet64')   # cervix, PVW


def hex_jac(P):
    d1 = (P[1] + P[2] + P[6] + P[5]) - (P[0] + P[3] + P[7] + P[4])
    d2 = (P[3] + P[2] + P[6] + P[7]) - (P[0] + P[1] + P[5] + P[4])
    d3 = (P[4] + P[5] + P[6] + P[7]) - (P[0] + P[1] + P[2] + P[3])
    return np.dot(d1, np.cross(d2, d3))


def _elems(t, name):
    m = re.search(r'<Elements type="hex8" name="%s">(.*?)</Elements>' % re.escape(name), t, re.S)
    return [[int(v) for v in c.split(',')] for c in re.findall(r'<elem id="\d+">([^<]+)</elem>', m.group(1))]


def _ray_tri(o, d, T):
    V0, V1, V2 = T[:, 0], T[:, 1], T[:, 2]
    A, B = V1 - V0, V2 - V0
    p = np.cross(d, B)
    det = (A * p).sum(1)
    ok = np.abs(det) > 1e-12
    inv = np.where(ok, 1 / np.where(ok, det, 1), 0)
    s = o - V0
    u = (s * p).sum(1) * inv
    q = np.cross(s, A)
    v = (q @ d) * inv
    tt = (B * q).sum(1) * inv
    hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (tt > 1e-6)
    return tt[hit].min() if hit.any() else None


def _closest_tri(p, T):
    """closest point of point p on triangles T (n, 3, 3): (index, point, barycentric weights)"""
    best = (np.inf, None, None, None)
    for i, (a, b, c) in enumerate(T):
        ab, ac, ap = b - a, c - a, p - a
        d1, d2 = ab @ ap, ac @ ap
        if d1 <= 0 and d2 <= 0:
            q, w = a, (1, 0, 0)
        else:
            bp = p - b; d3, d4 = ab @ bp, ac @ bp
            cp = p - c; d5, d6 = ab @ cp, ac @ cp
            vc, vb, va = d1 * d4 - d3 * d2, d5 * d2 - d1 * d6, d3 * d6 - d5 * d4
            if d3 >= 0 and d4 <= d3:
                q, w = b, (0, 1, 0)
            elif d6 >= 0 and d5 <= d6:
                q, w = c, (0, 0, 1)
            elif vc <= 0 and d1 >= 0 and d3 <= 0:
                v = d1 / (d1 - d3); q, w = a + v * ab, (1 - v, v, 0)
            elif vb <= 0 and d2 >= 0 and d6 <= 0:
                v = d2 / (d2 - d6); q, w = a + v * ac, (1 - v, 0, v)
            elif va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
                v = (d4 - d3) / ((d4 - d3) + (d5 - d6)); q, w = b + v * (c - b), (0, 1 - v, v)
            else:
                den = 1 / (va + vb + vc); v, u = vb * den, vc * den
                q, w = a + ab * v + ac * u, (1 - v - u, v, u)
        dd = np.linalg.norm(p - q)
        if dd < best[0]:
            best = (dd, i, q, w)
    return best


def close(t, K=3, c0=None, nn=None, maxd=12.0, mode='axis'):  # mode: axis (v3) | closest (v4) | end (v5)
    """mode 'axis' (v3, 2026-10-07): each cap node straight down the canal axis to the surface below it (fallback: the median
    distance) and the last ring tied node-to-node with its rest offset. mode 'closest' (v4, 2026-10-09, the user: the apex
    still shows a gap): each cap node's last ring lands ON its closest point of the cervix / PVW boundary (always defined) and
    is tied to that facet point (barycentric weights, zero rest offset), so the end face sits on the surface."""
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {d: _elems(t, d) for d in list(CURVES) + list(CAPWALLS) + ['_PickedSet347']}
    H = lambda p: (p - c0) @ nn
    # cervix / PVW boundary faces (triangles) and nodes
    cnt = Counter()
    for d in CAPWALLS:
        for el in E[d]:
            for fc in HF:
                cnt[frozenset(el[i] for i in fc)] += 1
    tris, trin, wn = [], [], set()
    for d in CAPWALLS:
        for el in E[d]:
            for fc in HF:
                f = [el[i] for i in fc]
                if cnt[frozenset(f)] == 1:
                    tris += [(X[f[0]], X[f[1]], X[f[2]]), (X[f[0]], X[f[2]], X[f[3]])]
                    trin += [(f[0], f[1], f[2]), (f[0], f[2], f[3])]
                    wn.update(f)
    T = np.array(tris)
    WN = np.array(sorted(wn))
    WX = np.array([X[n] for n in WN])
    a = t.index('<constraint name="canal_wrap_ends"')
    b = t.index('</constraint>', a)
    tied = {int(v) for v in re.findall(r'<node id="(\d+)" bc="x">1</node>', t[a:b])}
    wallnodes = {n for d in list(CAPWALLS) + ['_PickedSet347'] for el in E[d] for n in el}
    ccnt, cown = Counter(), {}
    for d in CURVES:
        for el in E[d]:
            for fi, fc in enumerate(HF):
                k = frozenset(el[i] for i in fc)
                ccnt[k] += 1
                cown[k] = (d, el, fi)
    nid = max(X) + 1
    new_nodes, new_hex, cons, log = [], {d + '_apex': [] for d in CURVES}, [], []
    for sgn in (1, -1):
        caps = []
        for k, c in ccnt.items():
            if c != 1:
                continue
            d, el, fi = cown[k]
            Q = np.array([X[el[i]] for i in HF[fi]])
            if sgn * Q.mean(0)[0] <= 1.5 or H(Q.mean(0)) > -70:
                continue
            hc = np.array([X[n] for n in el]).mean(0)
            nv = np.cross(Q[2] - Q[0], Q[3] - Q[1])
            nv /= np.linalg.norm(nv)
            if nv @ (Q.mean(0) - hc) < 0:
                nv = -nv
            if nv @ nn < -0.6:
                caps.append((d, el, fi))
        capn = set(n for d, el, fi in caps for n in (el[i] for i in HF[fi]))
        prev = {}
        for d, el, fi in caps:
            for i, j in EDGES:
                for p_, q_ in ((el[i], el[j]), (el[j], el[i])):
                    if p_ in capn and q_ not in capn:
                        prev[p_] = q_
        assert len(capn) == 35 and len(prev) == 35, (sgn, len(capn), len(prev))
        ring = {}
        D = {}
        for n in capn:
            dv = X[n] - X[prev[n]]
            step = np.linalg.norm(dv)
            dv /= step
            hit = _ray_tri(X[n], dv, T)
            D[n] = (dv, hit if hit is not None and hit < maxd else None, step)
        reach = [v[1] for v in D.values() if v[1] is not None]
        # target (2026-10-07 v3): straight down the canal axis (-nn, apical), each cap node to the cervix / PVW surface
        # below it (ray along -nn; fallback: the side's median distance): translated rings keep the hexes valid
        dn = -nn
        dist = {}
        for n in capn:
            hh = _ray_tri(X[n], dn, T)
            dist[n] = hh if hh is not None and hh < maxd else None
        got = [v for v in dist.values() if v is not None]
        med = float(np.median(got)) if got else 3.0
        tgt = {n: X[n] + dn * (dist[n] if dist[n] is not None else med) for n in capn}
        facet = {}
        if mode == 'end':
            # v5 (2026-10-09, the user: "the cap needs to go all the way to the end of the cervix"): every cap node straight
            # down the canal axis to the plane of the cervix's apical end; the last ring tied to its closest facet point
            hend = min(H(X[n]) for d in CAPWALLS[:1] for el in E[d] for n in el)
            L, stopped = [], 0
            for n in capn:
                full = H(X[n]) - hend
                hh = _ray_tri(X[n], dn, T)           # the axis ray meets the cervix / PVW first: stop on the surface
                if hh is not None and hh < full:
                    full, stopped = hh, stopped + 1
                tgt[n] = X[n] + dn * full
                dd, i, q, w = _closest_tri(tgt[n], T)
                facet[n] = (trin[i], w)
                L.append(full)
            log.append(f'{"right" if sgn > 0 else "left"} to the cervix end (H {hend:.1f}): extension {min(L):.1f}..{max(L):.1f} mm, '
                       f'{stopped} of {len(capn)} nodes stopped earlier on the cervix / PVW surface')
        if mode == 'closest':
            gaps = []
            for n in capn:
                dd, i, q, w = _closest_tri(X[n], T)
                tgt[n] = q; facet[n] = (trin[i], w); gaps.append(dd)
            log.append(f'{"right" if sgn > 0 else "left"} closest-point targets: all {len(capn)} cap nodes onto the cervix / '
                       f'PVW surface (distance median {np.median(gaps):.1f} mm, range {min(gaps):.1f}..{max(gaps):.1f})')
        log.append(f'{"right" if sgn > 0 else "left"} axis rays: {len(got)} of {len(capn)} hit (median {med:.1f} mm, '
                   f'range {min(got) if got else 0:.1f}..{max(got) if got else 0:.1f})')
        rings = [{}]
        for k in range(1, K + 1):
            rk = {}
            for n in capn:
                rk[n] = nid + len(new_nodes)
                new_nodes.append(X[n] + (tgt[n] - X[n]) * k / K)
            rings.append(rk)
        pos = lambda m: X[m] if m in X else new_nodes[m - nid]
        bad = 0
        for d, el, fi in caps:
            capf = set(el[i] for i in HF[fi])
            back = {prev[c]: c for c in capf}
            for k in range(1, K + 1):
                lo = (lambda c: c) if k == 1 else (lambda c, k=k: rings[k - 1][c])
                he = [rings[k][m] if m in capf else lo(back[m]) for m in el]
                if hex_jac(np.array([pos(m) for m in he])) <= 0:
                    bad += 1
                new_hex[d + '_apex'].append(he)
        ends = [n for n in capn if n in tied or n in wallnodes]
        for k in range(1, K + 1):
            for n in (capn if k == K else ends):
                m = rings[k][n]
                if k == K and mode in ('closest', 'end'):
                    cons.append((m, [(int(a), float(c)) for a, c in zip(*facet[n]) if c > 1e-9]))
                    continue
                w = WN[np.argmin(np.linalg.norm(WX - pos(m), axis=1))]
                cons.append((m, [(int(w), 1.0)]))
        hits = sum(v[1] is not None for v in D.values())
        log.append(f'{"right" if sgn > 0 else "left"}: {len(caps)} cap faces, {len(capn)} nodes, {hits} rays hit the cervix / '
                   f'PVW within {maxd:g} mm (median {np.median(reach) if reach else float("nan"):.1f} mm), {K} rings, '
                   f'{len(caps) * K} hexes, {bad} with non-positive Jacobian, ties: {len(ends)} end columns + last ring')
    # write
    blk = f'\t\t<Nodes name="apex_close_nodes">\n' + ''.join(
        f'\t\t\t<node id="{nid + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n' for i, p in enumerate(new_nodes)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    eblk, dom, damp = '', '', ''
    for dn, hx in new_hex.items():
        eblk += f'\t\t<Elements type="hex8" name="{dn}">\n'
        for he in hx:
            eblk += f'\t\t\t<elem id="{eid}">{",".join(map(str, he))}</elem>\n'
            eid += 1
        eblk += '\t\t</Elements>\n'
        base = dn[:-len('_apex')]
        dom += f'\t\t<SolidDomain name="{dn}" mat="{CURVES[base]}" />\n'
        cblk = re.search(r'<body_load name="wall_damping_%s"[^>]*>\s*<C>([^<]+)</C>' % re.escape(base), t)
        if cblk:
            damp += (f'\t\t<body_load name="wall_damping_{dn}" type="mass damping" elem_set="{dn}">\n'
                     f'\t\t\t<C>{cblk.group(1)}</C>\n\t\t</body_load>\n')
    i = t.index('\t</Mesh>')
    t = t[:i] + eblk + t[i:]
    i = t.index('\t</MeshDomains>')
    t = t[:i] + dom + t[i:]
    i = t.index('\t</Loads>')
    t = t[:i] + damp + t[i:]
    cx = '\t\t<constraint name="apex_close_ties" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n\t\t\t<penalty>10</penalty>\n' \
         '\t\t\t<maxaug>10</maxaug>\n' + ''.join(
             f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{m}" bc="{c}">1</node>\n'
             + ''.join(f'\t\t\t\t<node id="{w}" bc="{c}">{-cw:.6g}</node>\n' for w, cw in ws)
             + '\t\t\t</linear_constraint>\n' for m, ws in cons for c in 'xyz') + '\t\t</constraint>\n'
    i = t.index('\t</Constraints>')
    t = t[:i] + cx + t[i:]
    msg = ('NOT IN SOURCE (the user: a closed tube at the apex): each curve extruded beyond its apical end cap to the cervix / '
           f'PVW (tools/apex_close.py; {len(new_nodes)} nodes, {sum(len(v) for v in new_hex.values())} hexes in '
           f'canal_wrap_hex_avw_apex / _pvw_apex, {len(cons)} node ties apex_close_ties): ' + '; '.join(log))
    return t, msg
