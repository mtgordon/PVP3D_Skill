"""Batch 135 (2026-10-01, the user): "go back to the pre-narrowed version of the vaginal wall. Then make the width of the
entire vaginal wall 6mm narrower on each side and then attach the 3mm rounded edges to it. Have the same contacts as you
just ran with the levator. Make the thickness of the vaginal solid parts 2.75mm and make the curved pieces that thickness
as well"; then: "keep the top and bottom surfaces in the current location making the channel between them larger but
keeping the curve decently the same".
From L91_newline_vwyeoh_rhoi0_seamspr001 (the springs line before the tube):
  1. wall_reshape.reshape: every wall (AVW, cervix, PVW) 2.75 mm thick, its outer surface kept and its lumen surface
     moved out (the channel between the AVW and PVW opens by ~5 mm), the whole wall narrowed 6 mm on each side about its
     midline; the perineal body follows the shared PVW nodes (fading over 8 mm); the strip where the cervix curves into
     the PVW's end (no column) is filled smoothly from its neighbours. The canal seam springs are removed (the walls no
     longer meet there; the curved edges join them).
  2. the curved edges (wrap2): at every AVW / cervix seam station, a U of wall tissue joining the AVW's side face (its
     column A0 lumen edge .. A2 top) round laterally to the PVW's (P0 lumen edge .. P3 bottom), centred midway between
     A0 and P0: its outer surface runs from the AVW's top round to the PVW's bottom (radius ~ as before), its inner
     surface from the AVW's lumen edge to the PVW's (radius ~2.5 mm), so it is 2.75 mm thick like the walls. Rings at the
     union of the two columns' layer heights; 6 segments round; the U's ends share the column nodes where they exist and
     follow the columns elsewhere (linear constraints, penalty 10, maxaug 10). The AVW half (segments 0-2) is
     Vagina_AVW (canal_wrap_hex_avw), the PVW half Vagina_PVW (canal_wrap_hex_pvw).
  3. as the springs base: SlidingElastic1 seg_up 2, the PM as a structure (batch 113 step 2).
  4. the contacts of L134_tube_r30s_pm_narrow_wallsLA: the source's PVW_LA plus walls_LA (every other outer wall face,
     both wrap halves included, against the LA shells).
  L135_tube_thin275_narrow6   built only (the user looks at it first)
usage: py -3.10 build_batch135.py"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch90 import seam_nodes, _bedges  # noqa: E402
from build_batch118 import fold  # noqa: E402
from build_batch122 import hex_jac_ok  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM  # noqa: E402
from build_batch132 import seg_up2  # noqa: E402
from build_batch131 import facets_xml  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402
import wall_reshape  # noqa: E402
import build_batch133 as b133  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
NAME = 'L135_tube_thin275_narrow6'
NSEG = 6
EXTEND_DISTAL = False   # batch 142: carry the curves to the last seam station at the perineal-body end
EXTENDED = {}
LABEL = 'NOT IN SOURCE (2.75 mm walls narrowed 6 mm, curved edges)'


def wrap2(model, orig):
    """the curved edges on a reshaped model (model.src = the reshaped file); orig = the unreshaped file (for the seam
    pairing, which is topological)"""
    path = model.src
    f0 = Feb(path)
    seam, fP, _ = seam_nodes(Model7(orig))                # topological: taken from the unreshaped file
    X = {n: np.array(p, float) for n, p in f0.nodes.items()}
    cx, peb = set(f0.domain_nodes(tube.CX)), set(f0.domain_nodes(tube.PEB))
    rep, _ = fold(Model7(orig))                          # seam node <- PVW lumen-edge node, paired at rest
    inv = {a: p for p, a in rep.items()}
    pedge = sorted(inv.values())
    # exact through-thickness columns from the hex topology (wall_reshape.columns: lumen face -> opposite faces);
    # tube.columns' geometric walk can run off along the rim (it did at 4 stations a side here: the bumps)
    t0 = open(path, encoding='utf-8').read()
    lum = wall_reshape._lumen_faces(t0, 'SlidingElastic1Secondary') + wall_reshape._lumen_faces(t0, 'SlidingElastic1Primary')
    colA, colP = {}, {}
    for dm, dst in ((tube.AVW, colA), (tube.CX, colA), (tube.PVW, colP)):
        dn = set(f0.domain_nodes(dm))
        for n, col in wall_reshape.columns(f0, dm, [fc for fc in lum if all(q in dn for q in fc)]).items():
            if n not in dst or len(col) > len(dst[n]):
                dst[n] = col
    cons_ = model.root.find('Constraints')
    for c_ in list(cons_):
        if c_.get('name') == 'canal_seam_springs':
            cons_.remove(c_)
    new_nodes, cons, hexes, hexj = [], [], [], []
    nid0 = model.max_node_id() + 1

    def add(p):
        new_nodes.append(p)
        return nid0 + len(new_nodes) - 1

    def trim(col, n):
        return col[:n] if len(col) >= n else None
    report, rads = {}, []
    for sgn in (-1, 1):
        S = [a for a in seam if np.sign(X[a][0] - 1.5) == sgn]
        ax = np.linalg.svd(np.array([X[a] for a in S]) - np.mean([X[a] for a in S], axis=0))[2][0]
        S.sort(key=lambda a: X[a] @ ax)
        if EXTEND_DISTAL:
            # the end nearer the perineal body (lower z); its stations beyond the last paired one get partners from the
            # PVW lumen-edge nodes further along (PeB-shared ones allowed: the curves no longer merge nodes)
            fo = Feb(orig)
            Xo = {n: np.array(q, float) for n, q in fo.nodes.items()}
            cxo = set(fo.domain_nodes(tube.CX))
            eP = sorted({v for e in _bedges(fP) for v in e})
            idx = [i for i, a in enumerate(S) if a in inv]
            distal_end = len(S) - 1 if Xo[S[-1]][2] < Xo[S[0]][2] else 0
            ext = list(range(idx[-1] + 1, len(S))) if distal_end else list(range(idx[0] - 1, -1, -1))
            if ext:
                last_p = inv[S[idx[-1]] if distal_end else S[idx[0]]]
                axo = (Xo[S[ext[-1]]] - Xo[S[idx[-1] if distal_end else idx[0]]])
                axo /= np.linalg.norm(axo)
                used = set(inv.values())
                cand = [v for v in eP if v not in cxo and v not in used and np.sign(Xo[v][0] - 1.5) == sgn
                        and (Xo[v] - Xo[last_p]) @ axo > 0.2]
                cand.sort(key=lambda v: (Xo[v] - Xo[last_p]) @ axo)
                if cand:
                    # the end station takes its nearest free PVW lumen-edge node; stations in between are interpolated
                    end_a = S[ext[-1]]
                    p_ = min(cand, key=lambda v: np.linalg.norm(Xo[v] - Xo[end_a]))
                    if np.linalg.norm(Xo[p_] - Xo[end_a]) < 3.0:
                        inv[end_a] = p_
                        EXTENDED[end_a] = p_
        merged = [i for i, a in enumerate(S) if a in inv]
        if len(merged) < 2:
            continue
        rings = []
        for i in range(merged[0], merged[-1] + 1):
            a0 = S[i]
            A = trim(colA.get(a0, []), 3)
            if A is None:
                rings.append(None)
                continue
            if a0 in inv and inv[a0] in colP and trim(colP[inv[a0]], 4):
                Pn = trim(colP[inv[a0]], 4)
                Pp = [X[v] for v in Pn]
                pv = ('own', Pn)
            else:
                lo = max(j for j in merged if j < i)
                hi = min(j for j in merged if j > i)
                Pl, Ph = trim(colP.get(inv[S[lo]], []), 4), trim(colP.get(inv[S[hi]], []), 4)
                if Pl is None or Ph is None:
                    rings.append(None)
                    continue
                s = float(np.clip((X[a0] - X[S[lo]]) @ ax / max(1e-9, (X[S[hi]] - X[S[lo]]) @ ax), 0, 1))
                Pp = [(1 - s) * X[u] + s * X[w] for u, w in zip(Pl, Ph)]
                pv = ('between', Pl, Ph, s)
            Ap = [X[v] for v in A]
            c = 0.5 * (Ap[0] + Pp[0])                       # midway between the two lumen edges
            up, dn = Ap[-1] - c, Pp[-1] - c
            RA, RP = np.linalg.norm(up), np.linalg.norm(dn)
            up, dn = up / RA, dn / RP
            rA, rP = np.linalg.norm(Ap[0] - c), np.linalg.norm(Pp[0] - c)
            rads.append((rA, rP, RA, RP))
            sA = [(np.linalg.norm(p - c) - rA) / (RA - rA) for p in Ap]          # 0 .. 1 through the band
            sP = [(np.linalg.norm(p - c) - rP) / (RP - rP) for p in Pp]
            fr = sorted([(0.0, 'I'), (sA[1], 'A1'), (sP[1], 'P1'), (sP[2], 'P2'), (1.0, 'O')])
            t_ = ax - (ax @ up) * up
            w = np.cross(up, t_)
            w -= (w @ up) * up
            w /= np.linalg.norm(w)
            if w[0] * sgn < 0:
                w = -w

            def end_node(sv, col_nodes, col_s, kind_map):
                """a ring node on a column at band fraction sv: shared if a column node sits there, else constrained"""
                for k, ss in enumerate(col_s):
                    if abs(ss - sv) < 1e-9:
                        return col_nodes[k], None
                k = int(np.clip(np.searchsorted(col_s, sv) - 1, 0, len(col_s) - 2))
                u = (sv - col_s[k]) / (col_s[k + 1] - col_s[k])
                return None, [(col_nodes[k], 1 - u), (col_nodes[k + 1], u)]
            ring_ids = []
            for sv, kind in fr:
                seg_ids = []
                for j in range(NSEG + 1):
                    th = np.pi / 2 - np.pi * j / NSEG
                    if j == 0:
                        n_, wts = end_node(sv, A, sA, None)
                        if n_ is None:
                            p = Ap[0] + (Ap[-1] - Ap[0]) * 0  # placeholder, set below
                            k = int(np.clip(np.searchsorted(sA, sv) - 1, 0, 1))
                            u = (sv - sA[k]) / (sA[k + 1] - sA[k])
                            n_ = add((1 - u) * Ap[k] + u * Ap[k + 1])
                            cons.append((n_, wts))
                        seg_ids.append(n_)
                        continue
                    if j == NSEG:
                        if pv[0] == 'own':
                            n_, wts = end_node(sv, Pn, sP, None)
                            if n_ is None:
                                k = int(np.clip(np.searchsorted(sP, sv) - 1, 0, 2))
                                u = (sv - sP[k]) / (sP[k + 1] - sP[k])
                                n_ = add((1 - u) * Pp[k] + u * Pp[k + 1])
                                cons.append((n_, wts))
                            seg_ids.append(n_)
                            continue
                        _, Pl, Ph, s = pv
                        k = int(np.clip(np.searchsorted(sP, sv) - 1, 0, 2))
                        u = float(np.clip((sv - sP[k]) / (sP[k + 1] - sP[k]), 0, 1))
                        n_ = add((1 - u) * Pp[k] + u * Pp[k + 1])
                        cons.append((n_, [(Pl[k], (1 - s) * (1 - u)), (Pl[k + 1], (1 - s) * u),
                                          (Ph[k], s * (1 - u)), (Ph[k + 1], s * u)]))
                        seg_ids.append(n_)
                        continue
                    vert = up if th >= 0 else dn
                    d = np.cos(th) * w + abs(np.sin(th)) * vert
                    d /= np.linalg.norm(d)
                    tj = j / NSEG
                    r = (1 - tj) * (rA + sv * (RA - rA)) + tj * (rP + sv * (RP - rP))
                    seg_ids.append(add(c + r * d))
                ring_ids.append(seg_ids)
            rings.append(ring_ids)
        nh = 0
        P = lambda ids: np.array([X[i] if i in X else new_nodes[i - nid0] for i in ids])
        for i in range(len(rings) - 1):
            if rings[i] is None or rings[i + 1] is None:
                continue
            R0, R1 = rings[i], rings[i + 1]
            for j in range(NSEG):
                for k in range(len(R0) - 1):
                    hx = [R0[k][j], R0[k][j + 1], R0[k + 1][j + 1], R0[k + 1][j],
                          R1[k][j], R1[k][j + 1], R1[k + 1][j + 1], R1[k + 1][j]]
                    if hex_jac_ok(P(hx)) < 0:
                        hx = hx[4:] + hx[:4]
                    hexes.append(hx)
                    hexj.append(j)
                    nh += 1
        report[sgn] = (sum(r is not None for r in rings), len(rings), nh)
    ids = model.add_nodes('canal_wrap_nodes', new_nodes)
    assert ids[0] == nid0
    blocks = model.mesh.findall('Elements')
    start = max(int(e.get('id')) for b in blocks for e in b) + 1
    i0 = list(model.mesh).index(blocks[-1])
    parts = [('canal_wrap_hex_avw', 'Vagina_AVW', [h for h, j in zip(hexes, hexj) if j < NSEG / 2]),
             ('canal_wrap_hex_pvw', 'Vagina_PVW', [h for h, j in zip(hexes, hexj) if j >= NSEG / 2])]
    q0 = 0
    doms = model.root.find('MeshDomains')
    for q, (nm, mat, hs) in enumerate(parts):
        eh = ET.Element('Elements', type='hex8', name=nm)
        for qq, h in enumerate(hs):
            ET.SubElement(eh, 'elem', id=str(start + q0 + qq)).text = ','.join(map(str, h))
        q0 += len(hs)
        model.mesh.insert(i0 + 1 + q, eh)
        ET.SubElement(doms, 'SolidDomain', name=nm, mat=mat)
    c = ET.SubElement(cons_, 'constraint', name='canal_wrap_ends', type='linear constraint')
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
    allP = {**X, **{nid0 + i: p for i, p in enumerate(new_nodes)}}
    jh = np.array([hex_jac_ok(np.array([allP[n] for n in h])) for h in hexes])
    rads = np.array(rads)
    model.log.append(f'{LABEL}: curved edges: per side (stations wrapped / spanned, hexes) {report}; {len(new_nodes)} new '
                     f'nodes, {len(hexes)} hex8 (canal_wrap_hex_avw Vagina_AVW / canal_wrap_hex_pvw Vagina_PVW); inner '
                     f'radius AVW / PVW end median {np.median(rads[:, 0]):.2f} / {np.median(rads[:, 1]):.2f} mm, outer '
                     f'{np.median(rads[:, 2]):.2f} / {np.median(rads[:, 3]):.2f} mm; {len(cons)} end nodes following the '
                     f'columns (canal_wrap_ends); hex Jacobians all positive: {bool((jh > 0).all())}')
    return jh


def walls_la(model):
    """the contacts of L134: walls_LA (every outer wall face not in PVW_LA_primary, both wrap halves included) against
    the LA (PVW_LA_secondary), settings as PVW_LA"""
    tmp = os.path.join(SCRATCH, NAME + '_pre_contact.feb')
    model.write(tmp)
    f, t = Feb(tmp), open(tmp, encoding='utf-8').read()

    def surf(name):
        i = t.index(f'<Surface name="{name}">')
        blk = t[i:t.index('</Surface>', i)]
        return [tuple(int(v) for v in m.split(',')) for m in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', blk)]
    WALLS = (tube.AVW, tube.CX, tube.PVW, 'canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
    lumen = {tuple(sorted(fc)) for s in ('SlidingElastic1Primary', 'SlidingElastic1Secondary') for fc in surf(s)}
    pvwla = {tuple(sorted(fc)) for fc in surf('PVW_LA_primary')}
    wn = {n for dm in WALLS for n in f.domain_nodes(dm)}
    faces = [fc for fc in tube.boundary_faces(f, WALLS + (tube.PEB,))
             if all(n in wn for n in fc) and tuple(sorted(fc)) not in lumen and tuple(sorted(fc)) not in pvwla]
    inside = [n for n in {n for fc in faces for n in fc}
              if (lambda r: r[0] < b133.FAR and r[1] < 2.0)(b133.signed_la(np.array(f.nodes[n])))]
    dom_of = {}
    for dm in WALLS:
        for n in f.domain_nodes(dm):
            dom_of.setdefault(n, dm)
    srf = ET.fromstring(facets_xml('walls_LA_primary', faces))
    mesh = model.mesh
    anchor = next(x for x in mesh if x.tag == 'SurfacePair' and x.get('name') == 'PVW_LA')   # after the LA surfaces
    mesh.insert(list(mesh).index(anchor), srf)
    sp = ET.Element('SurfacePair', name='walls_LA')
    ET.SubElement(sp, 'primary').text = 'walls_LA_primary'
    ET.SubElement(sp, 'secondary').text = 'PVW_LA_secondary'
    mesh.insert(list(mesh).index(anchor), sp)
    tmpl = next(x for x in model.root.find('Contact') if x.get('name') == 'PVW_LA')
    import copy
    con = copy.deepcopy(tmpl)
    con.set('name', 'walls_LA')
    con.set('surface_pair', 'walls_LA')
    model.root.find('Contact').append(con)
    model.log.append(f'NOT IN SOURCE: contact walls_LA, {len(faces)} outer wall faces not in PVW_LA_primary '
                     f'({dict(Counter(dom_of[fc[0]] for fc in faces))}) against the LA shells, as PVW_LA; wall nodes '
                     f'starting inside the LA\'s thickness: {len(inside)}')
    return inside


if __name__ == '__main__':
    os.makedirs(SCRATCH, exist_ok=True)
    assert not os.path.exists(os.path.join(RUNS, NAME)), f'{NAME} exists; not overwriting'
    orig = os.path.join(RUNS, BASE, BASE + '.feb')
    m0 = Model7(orig)
    rep = wall_reshape.reshape(m0)
    print('reshape:', rep)
    mid = os.path.join(SCRATCH, NAME + '_reshaped.feb')
    m0.write(mid)
    m = Model113(mid)
    m.log.append(f'{LABEL}: wall_reshape: {rep}')
    jh = wrap2(m, orig)
    print(m.log[-1])
    seg_up2(m)
    step2(INNER + OUTER_BOTTOM)(m)
    inside = walls_la(m)
    print(m.log[-1])
    emit(NAME, m)
