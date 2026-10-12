"""Extend the canal walls past the seam (2026-10-07, the user: "extend the avw, pvw, and cervix out 3mm but do not apply a
load to any of those elements"). On the contact-only seam (SlidingElastic1: AVW + cervix secondary, PVW + PeB primary) the
soft walls slide until an edge node runs off the other wall's last facet and flips on and off it (the stall; up to 5.4 mm
past the edge in L167_pm_outer_band3x3_damp_P2 before it stalled, tools/seam_slide.py). An unloaded lip on each wall gives
the edge nodes facets to slide onto, without stiffening the seam.
How: each wall's lateral seam side (the boundary hex faces through the wall thickness under each along-canal seam edge of
the lumen surface; AVW 2 layers, PVW 3) is extruded outward in K rings to L mm. Every node of a column moves by the same
vector, so the wall thickness is kept; the direction is shared by the two walls at each point of the seam (the average of
the AVW / cervix side normal and the facing PVW / PeB side normal, with the along-canal part removed), so the two lips stay
stacked at the seam's rest gap instead of crossing. The new hexes go in domains <wall>_ext with the wall's material and
damping. Their lumen faces join the seam contact surfaces (SlidingElastic1Secondary / Primary) only: no pressure, no other
contact. NOT IN SOURCE.
usage (module): t2, msg = extend(t, L=3.0, K=2)"""
import re
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]   # outward hex faces
WALLS = {'_PickedSet347': 'AVW', '_PickedSet346': 'cervix', '_PickedSet64': 'PVW', '_PickedSet66': 'PeB'}
SIDE = {'AVW': 'SlidingElastic1Secondary', 'cervix': 'SlidingElastic1Secondary',
        'PVW': 'SlidingElastic1Primary', 'PeB': 'SlidingElastic1Primary'}


def _unit(v):
    return v / np.linalg.norm(v)


def _qn(P):
    """quad normal (right-hand rule over the node order)"""
    return _unit(np.cross(P[2] - P[0], P[3] - P[1]))


def _bedges(faces):
    c = defaultdict(int)
    for f in faces:
        for k in range(len(f)):
            c[tuple(sorted((f[k], f[(k + 1) % len(f)])))] += 1
    return [e for e, n in c.items() if n == 1]


def _dseg(p, a, b):
    t = np.clip((p - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1)
    return np.linalg.norm(p - a - t * (b - a))


def hex_jac_min(P):
    """smallest corner Jacobian determinant of a hex8"""
    nb = [(1, 3, 4), (2, 0, 5), (3, 1, 6), (0, 2, 7), (5, 7, 0), (6, 4, 1), (7, 5, 2), (4, 6, 3)]
    out = []
    for i, (a, b, c) in enumerate(nb):
        J = np.array([P[a] - P[i], P[b] - P[i], P[c] - P[i]])
        out.append(np.linalg.det(J) * (1 if i < 4 else -1))
    return min(out)


def extend(t, L=3.0, K=2, tol=0.75):
    root = ET.fromstring(t)
    mesh = root.find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')
         if e.get('name') in WALLS}
    surf = {s.get('name'): [tuple(int(v) for v in f.text.split(',')) for f in s] for s in mesh.findall('Surface')}
    fS, fP = surf['SlidingElastic1Secondary'], surf['SlidingElastic1Primary']
    # every wall hex face, and which faces are on the boundary of all the walls together
    fcount = defaultdict(int)
    for w, els in E.items():
        for el in els:
            for hf in HF:
                fcount[frozenset(el[i] for i in hf)] += 1
    # seam edges: lumen boundary edges whose nodes lie within tol of the other surface's boundary edge, along the canal
    C = np.array([X[v] for f in fS for v in f])
    ax = np.linalg.svd(C - C.mean(0), full_matrices=False)[2][0]
    bS, bP = _bedges(fS), _bedges(fP)
    def seam(own, other):
        near = lambda v: min(_dseg(X[v], X[a], X[b]) for a, b in other) < tol  # noqa: E731
        ok = {v for e in own for v in e}
        ok = {v for v in ok if near(v)}
        return [e for e in own if e[0] in ok and e[1] in ok and abs(_unit(X[e[1]] - X[e[0]]) @ ax) > 0.5]
    sS, sP = seam(bS, bP), seam(bP, bS)
    # each seam edge -> its wall, the lumen hex, and the column of side faces through the thickness
    lumen = set(map(frozenset, fS)) | set(map(frozenset, fP))
    owner = {}
    for w, els in E.items():
        for el in els:
            for hf in HF:
                fs = frozenset(el[i] for i in hf)
                if fs in lumen:
                    owner.setdefault(fs, []).append((w, el))
    bfaces = defaultdict(list)          # boundary wall face (ordered outward) by edge
    for w, els in E.items():
        for el in els:
            for hf in HF:
                q = tuple(el[i] for i in hf)
                if fcount[frozenset(q)] == 1:
                    for k in range(4):
                        bfaces[frozenset((q[k], q[(k + 1) % 4]))].append((w, q))
    cols, log = [], []
    for side, edges, faces in (('sec', sS, fS), ('pri', sP, fP)):
        for a, b in edges:
            lf = [f for f in faces if a in f and b in f][0]
            w, el = owner[frozenset(lf)][0]
            # the first side face: a boundary face of the same hex holding a, b (not the lumen face)
            cand = [q for ww, q in bfaces[frozenset((a, b))] if set(q) <= set(el) and frozenset(q) != frozenset(lf)]
            if not cand:
                log.append(f'{a}-{b}: no side face')
                continue
            q = cand[0]
            col = [q]
            top = {a, b}
            for _ in range(4):
                far = [v for v in q if v not in top]
                nq = _qn(np.array([X[v] for v in q]))
                nxt = [qq for ww, qq in bfaces[frozenset(far)] if frozenset(qq) != frozenset(q)
                       and _qn(np.array([X[v] for v in qq])) @ nq > 0.5]
                if not nxt:
                    break
                top, q = set(far), nxt[0]
                col.append(q)
            cols.append((side, WALLS[w], w, (a, b), col))
    # the outward direction per lumen seam node: own side normal + the facing wall's, along-canal part removed
    nrm = defaultdict(list)
    for side, wn, w, (a, b), col in cols:
        n = _qn(np.array([X[v] for v in col[0]]))
        nrm[(side, a)].append(n)
        nrm[(side, b)].append(n)
    own_n = {k: _unit(np.mean(v, axis=0)) for k, v in nrm.items()}
    secn = [v for (s, v) in own_n if s == 'sec']
    prin = [v for (s, v) in own_n if s == 'pri']
    D = {}
    for (s, v), n in own_n.items():
        others = prin if s == 'sec' else secn
        o = min(others, key=lambda u: np.linalg.norm(X[u] - X[v]))
        d = n + own_n[('pri' if s == 'sec' else 'sec', o)]
        d -= (d @ ax) * ax
        D[v] = _unit(d)
    # every column node follows its lumen node's direction (the nearest lumen seam node of the column)
    nid = max(X) + 1
    newn = {}                            # (node, ring) -> id
    coords = []
    hexes = defaultdict(list)
    lumen_new = defaultdict(list)        # surface -> new faces
    shear = []
    def node(v, k, d):
        nonlocal nid
        if k == 0:
            return v
        if (v, k) not in newn:
            newn[(v, k)] = nid
            coords.append(X[v] + d * L * k / K)
            X[nid] = coords[-1]
            nid += 1
        return newn[(v, k)]
    colnode_dir = {}
    for side, wn, w, (a, b), col in cols:
        for q in col:
            for v in q:
                lv = min((a, b), key=lambda u: np.linalg.norm(X[u] - X[v]))
                colnode_dir.setdefault(v, D[lv])
    bad = 0
    for side, wn, w, (a, b), col in cols:
        for qi, q in enumerate(col):
            shear.append(np.degrees(np.arccos(np.clip(_qn(np.array([X[v] for v in q])) @ colnode_dir[q[0]], -1, 1))))
            for k in range(K):
                lo = [node(v, k, colnode_dir[v]) for v in q]
                hi = [node(v, k + 1, colnode_dir[v]) for v in q]
                h = lo + hi
                if hex_jac_min(np.array([X[v] for v in h])) <= 0:
                    h = lo[::-1] + hi[::-1]
                    if hex_jac_min(np.array([X[v] for v in h])) <= 0:
                        bad += 1
                hexes[w].append(h)
                if qi == 0:   # the lip's lumen face: the new hex face holding a, b (ring k) and their ring k + 1 copies
                    want = {node(a, k, colnode_dir[a]), node(b, k, colnode_dir[b]),
                            node(a, k + 1, colnode_dir[a]), node(b, k + 1, colnode_dir[b])}
                    f = [tuple(h[i] for i in hf) for hf in HF if {h[i] for i in hf} == want]
                    lumen_new[SIDE[wn]].append(f[0])
    # the lips' lumen faces must face the same way as the seam's own lumen faces
    for sn, fl in lumen_new.items():
        fixed = []
        for f in fl:
            # orient like the nearest existing lumen face
            P = np.array([X[v] for v in f])
            c = P.mean(0)
            nf = min(surf[sn], key=lambda g: np.linalg.norm(np.mean([X[v] for v in g], axis=0) - c))
            fixed.append(f if _qn(P) @ _qn(np.array([X[v] for v in nf])) > 0 else f[::-1])
        lumen_new[sn] = fixed
    # gap between the two lips at the outer ring (should stay at the seam's rest gap)
    outer = {v: X[newn[(v, K)]] for (v, k) in newn if k == K and v in D}
    so = [outer[v] for v in secn if v in outer]
    po = [outer[v] for v in prin if v in outer]
    gaps = [min(np.linalg.norm(p - s) for p in po) for s in so] if so and po else [float('nan')]
    # write
    blk = '\t\t<Nodes name="seam_ext_nodes">\n' + ''.join(
        f'\t\t\t<node id="{max(X) - len(coords) + 1 + i}">{p[0]:.7g},{p[1]:.7g},{p[2]:.7g}</node>\n'
        for i, p in enumerate(coords)) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    mats = {d.get('name'): d.get('mat') for d in root.find('MeshDomains')}
    eblk, dom, damp, cnt = '', '', '', []
    for w, hx in hexes.items():
        dn = f'{WALLS[w]}_ext'
        eblk += f'\t\t<Elements type="hex8" name="{dn}">\n'
        for h in hx:
            eblk += f'\t\t\t<elem id="{eid}">{",".join(map(str, h))}</elem>\n'
            eid += 1
        eblk += '\t\t</Elements>\n'
        dom += f'\t\t<SolidDomain name="{dn}" mat="{mats[w]}" />\n'
        cb = re.search(r'<body_load name="wall_damping_%s"[^>]*>\s*<C>([^<]+)</C>' % re.escape(w.lstrip('_')), t)
        if cb:
            damp += (f'\t\t<body_load name="wall_damping_{dn}" type="mass damping" elem_set="{dn}">\n'
                     f'\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n')
        cnt.append(f'{dn} {len(hx)} ({mats[w]}{", damped" if cb else ""})')
    i = t.index('\t</Mesh>')
    t = t[:i] + eblk + t[i:]
    i = t.index('\t</MeshDomains>')
    t = t[:i] + dom + t[i:]
    i = t.index('\t</Loads>')
    t = t[:i] + damp + t[i:]
    for sn, fl in lumen_new.items():
        a = t.index(f'<Surface name="{sn}">')
        b = t.index('</Surface>', a)
        fid = max(int(v) for v in re.findall(r'<quad4 id="(\d+)"', t[a:b])) + 1
        add = ''.join(f'\t\t\t<quad4 id="{fid + j}">{",".join(map(str, f))}</quad4>\n' for j, f in enumerate(fl))
        t = t[:b] + add.lstrip('\t') + '\t\t' + t[b:] if t[b - 1] == '\t' and False else t[:b].rstrip('\t') + add + '\t\t' + t[b:]
    msg = (f'NOT IN SOURCE (the user: extend the walls past the seam, unloaded): each wall\'s lateral seam side extruded {L:g} mm '
           f'outward in {K} rings (tools/seam_extend.py): {len(sS)} AVW / cervix and {len(sP)} PVW / PeB seam edges, '
           f'{len(coords)} nodes, {"; ".join(cnt)}; {bad} hexes with a non-positive Jacobian; side-face tilt from the shared '
           f'direction median {np.median(shear):.0f} deg, max {np.max(shear):.0f}; lip-to-lip gap at the outer ring median '
           f'{np.median(gaps):.2f} mm, max {np.max(gaps):.2f}; lumen faces added to the seam contact only '
           f'({", ".join(f"{k} +{len(v)}" for k, v in lumen_new.items())}); no pressure on them'
           + (f'; skipped: {"; ".join(log)}' if log else ''))
    return t, msg


if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from variants7 import JOBS
    from build_batch162 import write
    RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
    src, name, L, K = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    t2, msg = extend(t, L, K)
    print(msg)
    if '--dry' not in sys.argv:
        notes = open(os.path.join(RUNS, src, src + '.feb.changes.txt'), encoding='utf-8').read()
        write(name, t2, notes + msg + f'\n(base {src})\n')
