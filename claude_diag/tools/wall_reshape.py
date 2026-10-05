"""Reshape the vaginal walls (2026-10-01, the user, for batch 135): narrow the whole wall by NARROW mm on each side and make
every wall solid THICK mm thick, keeping the outer surfaces (the AVW's top, the PVW's bottom, the cervix's outside) where
they are, so the channel between the walls gets larger.
  * columns: each wall (AVW, cervix, PVW; hex8, layered through the thickness: 3 / 3 / 4 nodes) is walked from its lumen
    faces (the faces of SlidingElastic1Secondary / Primary on that domain) to the opposite face of each hex, layer by
    layer, giving every lumen node's column of nodes to the outer surface;
  * thinning: in each column the outer node stays, the others move toward it so the column is THICK long
    (c_k' = c_out + (c_k - c_out) THICK / T);
  * narrowing: every wall node's lateral offset from the midline x0 (the seam's mean x) is scaled by (h - NARROW) / h,
    h = the lateral half-width at that point (|x - x0| of the nearest seam node, nearest in y-z), so each side's edge
    comes in NARROW mm and the midline stays;
  * the perineal body shares the PVW's distal end face: its nodes follow the shared nodes' displacement, fading
    linearly to zero FADE mm away from them, so its elements are not torn.
  * flare (optional, batch 136): {side sign: (junction node ids, L)}: the PVW's narrowing fades out toward the perineal body,
    from full where the curved edges stop (L mm from the PVW / PeB junction) to none at the junction (smoothstep in the
    distance to the junction nodes), and the perineal body follows only the thinning, not the narrowing.
usage: from wall_reshape import reshape; reshape(model) on a Model7 (edits node coordinates in place; returns a report)"""
import re

import numpy as np

import tube
from build_batch90 import seam_nodes
from variants7 import Model7

NARROW, THICK, FADE = 6.0, 2.75, 8.0          # batch 137 sets NARROW = 3.0 before reshaping
HEXF = tube.HEXF
OPP = {0: 1, 1: 0, 2: 4, 3: 5, 4: 2, 5: 3}       # opposite faces of tube.HEXF (bottom/top, and the side pairs)


def _lumen_faces(t, name):
    i = t.index(f'<Surface name="{name}">')
    blk = t[i:t.index('</Surface>', i)]
    return [tuple(int(v) for v in m.split(',')) for m in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', blk)]


def columns(f, dom, lumen):
    """inner node -> [nodes from the inner face to the outer surface] for a layered hex domain: from each inner face,
    step to the opposite face of its hex, layer by layer"""
    els = f.elem_blocks[dom][1]
    face_of = {}
    for e, c in els.items():
        for k, fc in enumerate(HEXF):
            face_of.setdefault(tuple(sorted(c[i] for i in fc)), []).append((e, k))
    cols = {}
    for lf in lumen:
        cur = tuple(sorted(lf))
        if cur not in face_of:
            continue
        seq, seen = [], set()
        while cur in face_of:
            nxt = [(e, k) for e, k in face_of[cur] if e not in seen]
            if not nxt:
                break
            e, k = nxt[0]
            seen.add(e)
            c = els[e]
            a, b = [c[i] for i in HEXF[k]], [c[i] for i in HEXF[OPP[k]]]
            pairs = {}
            for i, j in tube.HEXE:
                if c[i] in a and c[j] in b:
                    pairs[c[i]] = c[j]
                elif c[j] in a and c[i] in b:
                    pairs[c[j]] = c[i]
            seq.append(pairs)
            cur = tuple(sorted(b))
        for n in lf:
            col = [n]
            for pairs in seq:
                if col[-1] in pairs:
                    col.append(pairs[col[-1]])
            if n not in cols or len(col) > len(cols[n]):
                cols[n] = col
    return cols


def grow_inner(fb, dom, seeds, exclude_nodes, max_angle=45.0):
    """region-grow the inner boundary sheet of a layered domain from seed faces over boundary faces sharing an edge and
    turning less than max_angle; faces made only of exclude_nodes (interfaces with other walls) are skipped"""
    bf = tube.boundary_faces(fb, (dom,))
    key = lambda fc: tuple(sorted(fc))
    X = fb.nodes
    def nrm(fc):
        P = np.array([X[n] for n in fc])
        v = np.cross(P[2] - P[0], P[3] - P[1]) if len(fc) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
        return v / np.linalg.norm(v)
    edges = {}
    for fc in bf:
        for a, b in zip(fc, fc[1:] + fc[:1]):
            edges.setdefault(tuple(sorted((a, b))), []).append(key(fc))
    byk = {key(fc): fc for fc in bf}
    inner = {key(fc) for fc in seeds if key(fc) in byk}
    stack = list(inner)
    while stack:
        k = stack.pop()
        fc = byk[k]
        n0 = nrm(fc)
        for a, b in zip(fc, fc[1:] + fc[:1]):
            for k2 in edges[tuple(sorted((a, b)))]:
                if k2 in inner or all(n in exclude_nodes for n in k2):
                    continue
                if abs(n0 @ nrm(byk[k2])) >= np.cos(np.radians(max_angle)):
                    inner.add(k2)
                    stack.append(k2)
    return [byk[k] for k in inner]


def reshape(model, flare=None, flare_doms=None):
    f = model                                       # Model7 offers nodes(); we need the Feb view for topology
    from febmodel import Feb
    fb = Feb(model.src)
    t = open(model.src, encoding='utf-8').read()
    X = {n: np.array(p, float) for n, p in fb.nodes.items()}
    sec, pri = _lumen_faces(t, 'SlidingElastic1Secondary'), _lumen_faces(t, 'SlidingElastic1Primary')
    walls = (tube.AVW, tube.CX, tube.PVW)
    targets = {}
    report = {}
    for dom in walls:
        dn = set(fb.domain_nodes(dom))
        lum = [fc for fc in sec + pri if all(n in dn for n in fc)]
        cols = columns(fb, dom, lum)
        fixed_inner = {}
        if dom == tube.CX:
            # the cervix's inner face also continues past its lumen faces (no contact surface there): grow it across
            # neighbouring boundary faces that bend less than 45 degrees (the rim faces turn ~90 degrees)
            inner = grow_inner(fb, dom, lum, set(fb.domain_nodes(tube.PVW)) | set(fb.domain_nodes(tube.AVW)))
            cols = columns(fb, dom, inner)
        L = [len(c) for c in cols.values()]
        T = [np.linalg.norm(X[c[-1]] - X[c[0]]) for c in cols.values()]
        covered = {n for c in cols.values() for n in c}
        report[dom] = (len(cols), sorted(set(L)), round(float(np.median(T)), 2), len(dn - covered))
        for n0, c in cols.items():
            if n0 in fixed_inner and n0 not in {m for fc in lum for m in fc}:
                inn = X[c[0]]
                Tc = np.linalg.norm(X[c[-1]] - inn)
                for n in c[1:]:
                    targets.setdefault(n, []).append(inn + (X[n] - inn) * THICK / Tc)
                continue
            out = X[c[-1]]
            targets.setdefault(c[-1], []).append(out)          # the outer surface stays
            Tc = np.linalg.norm(X[c[0]] - out)
            for n in c[:-1]:
                targets.setdefault(n, []).append(out + (X[n] - out) * THICK / Tc)
    # thinning displacement (nodes reached from two domains: the mean)
    disp = {n: np.mean(v, axis=0) - X[n] for n, v in targets.items()}
    # wall nodes in no column (the strip where the cervix curves round into the PVW's end): their thinning displacement
    # is filled in smoothly from their neighbours (Jacobi averaging over the hex edges, the column nodes held)
    wall_all = {n for dm in walls for n in fb.domain_nodes(dm)}
    free = [n for n in wall_all if n not in disp]
    adj = tube.adjacency(fb, walls)
    cur = {n: np.zeros(3) for n in free}
    for _ in range(400):
        nxt = {}
        for n in free:
            nb = [disp[m] if m in disp else cur[m] for m in adj[n] if m in disp or m in cur]
            nxt[n] = np.mean(nb, axis=0) if nb else np.zeros(3)
        cur = nxt
    disp.update(cur)
    report['filled smoothly (no column)'] = len(free)
    # narrowing (applied to the thinned positions)
    seam, _, _ = seam_nodes(Model7(model.src))
    S = np.array([X[s] for s in seam])
    x0 = float(np.mean(S[:, 0]))
    wall_nodes = {n for dm in walls for n in fb.domain_nodes(dm)}
    pvw_nodes = set(fb.domain_nodes(tube.PVW))
    # the walls that flare (batch 136: the PVW only; batch 142: all three, so the curves on the walls flare with them)
    flare_nodes = {n for dm in (flare_doms or (tube.PVW,)) for n in fb.domain_nodes(dm)}
    if flare is not None:
        flare = {k: (np.array([X[j] for j in v[0]]), v[1]) for k, v in flare.items()}
    new = {}
    for n in wall_nodes:
        p = X[n] + disp.get(n, 0)
        k = int(np.argmin(np.linalg.norm(S[:, 1:] - X[n][1:], axis=1)))
        h = abs(S[k, 0] - x0)
        sc = max(h - NARROW, 0.25 * h) / h
        if flare is not None and n in flare_nodes:
            side = 1 if X[n][0] >= x0 else -1
            J, Lf = flare[side]
            u = float(np.clip(np.linalg.norm(J - X[n], axis=1).min() / Lf, 0, 1))
            wgt = u * u * (3 - 2 * u)                       # smoothstep: 0 at the junction, 1 at the wrap's end
            sc = 1 - wgt * (1 - sc)
        p = p.copy()
        p[0] = x0 + (p[0] - x0) * sc
        new[n] = p
    # the perineal body follows the shared nodes, fading out
    peb = set(fb.domain_nodes(tube.PEB))
    shared = [n for n in peb if n in new]
    SH = np.array([X[n] for n in shared])
    DSH = np.array([new[n] - X[n] for n in shared])
    moved_peb = 0
    for n in peb - set(shared):
        d = np.linalg.norm(SH - X[n], axis=1)
        if d.min() >= FADE:
            continue
        w = 1.0 / np.maximum(d, 1e-6) ** 2
        u = (w[:, None] * DSH).sum(0) / w.sum() * (1 - d.min() / FADE)
        new[n] = X[n] + u
        moved_peb += 1
    model.set_node_coords(new) if hasattr(model, 'set_node_coords') else _set_coords(model, new)
    return {'columns per wall (lumen nodes, column lengths, median thickness, wall nodes not in a column)': report,
            'x0': round(x0, 2), 'PeB nodes moved': moved_peb, 'wall nodes moved': len(wall_nodes)}


def _set_coords(model, new):
    for blk in model.mesh.findall('Nodes'):
        for nd in blk:
            n = int(nd.get('id'))
            if n in new:
                p = new[n]
                nd.text = f'{p[0]:.8g},{p[1]:.8g},{p[2]:.8g}'
