"""Batch 413 (2026-10-22; the user's pictures of the strap at the curve end: "it should fill in where I added the red here (and the
same on the other side)" -> proposal strap_wide_proposal.png -> "Let's have it go another 2.5mm (or so) along the pm arch").
NOT IN SOURCE (the strap is ours): each fan strap (PM_avw_bottom_left/right_sheet_solid, one layer of penta6) is widened on its
MEDIAL side by a strip (new domain <strap>_wide, the strap's material, fibre directions copied from the nearest strap element):
- its strap-side edge = the strap's medial boundary (curve corner -> PM-end anchor), sharing those bottom and top nodes;
- its curve edge = the curve's end-face rim from the strap's corner to the first curve node on the AVW (_PickedSet347): the rim
  nodes are used directly, the other edge nodes are tied onto their rim edge (linear constraints, strap_wide_ties), as the
  strap's own curve edge (direct tie);
- its PM-end edge = ARCH mm along the PM_Plane inner edge from the strap's last anchor: new nodes fixed with the strap's anchors
  (NodeSet curve_end_PM_anchors, zero displacement);
- the far edge straight from the AVW node to the arch end; interior by a Coons patch; NW divisions across, one row per strap
  medial edge segment; top layer = bottom + the strap's thickness vector at that row.
Surfaces: wherever the strap's bottom / top triangles are, the strip's are added; side quads of the strap on the medial edge
(now inside) are removed and the strip's far-edge side quads added where the strap has side quads.
usage: py -3.10 build_batch413.py SRC_RUN NEW_RUN [--arch 5.0] [--nw 5]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch278 import lc  # noqa: E402
from build_batch162 import write  # noqa: E402

STRAPS = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
HF8 = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
W3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def wvol(P):
    return sum(np.dot(np.cross(P[b] - P[a], P[d] - P[a]), P[e] - P[a]) / 6 for a, b, d, e in W3)


def resample(poly, n):
    """n + 1 points equally spaced by arc length along a polyline; also (segment index, fraction) of each"""
    poly = np.asarray(poly); seg = np.linalg.norm(np.diff(poly, axis=0), axis=1); cum = np.r_[0, np.cumsum(seg)]
    out = []
    for s in np.linspace(0, cum[-1], n + 1):
        i = min(np.searchsorted(cum, s, side='right') - 1, len(seg) - 1)
        f = (s - cum[i]) / seg[i] if seg[i] > 0 else 0.0
        out.append((poly[i] + f * (poly[i + 1] - poly[i]), i, f))
    return out


def main(src, new, arch=5.0, nw=5):
    t, X, E = load(src)
    X = {k: np.asarray(p, float) for k, p in dict(X).items()}
    dn = defaultdict(set)
    for k, (n, c) in E.items():
        for v in c:
            dn[v].add(n)
    iscurve = lambda v: any(d.startswith('canal_wrap') for d in dn[v])   # noqa: E731
    m = re.search(r'<NodeSet name="curve_end_PM_anchors">([^<]*)</NodeSet>', t)
    anchors = {int(v) for v in re.findall(r'\d+', m.group(1))}
    # curve end-face rim edges (boundary edges of the distal end faces)
    cur = {k: c for k, (n, c) in E.items() if n.startswith('canal_wrap') and len(c) == 8}
    cnt = Counter(frozenset(c[i] for i in f) for c in cur.values() for f in HF8)
    endf = [[c[i] for i in f] for c in cur.values() for f in HF8[:2]
            if cnt[frozenset(c[i] for i in f)] == 1 and np.mean([X[c[i]][2] for i in f]) < -37]
    ec = Counter(frozenset((f[i], f[(i + 1) % 4])) for f in endf for i in range(4))
    rim = defaultdict(set)
    for e, n in ec.items():
        if n == 1:
            a, b = tuple(e); rim[a].add(b); rim[b].add(a)
    # PM_Plane boundary loop
    pq = [c for n, c in E.values() if n == 'PM_Plane']
    pe = Counter(frozenset((q[i], q[(i + 1) % 4])) for q in pq for i in range(4))
    padj = defaultdict(set)
    for e, n in pe.items():
        if n == 1:
            a, b = tuple(e); padj[a].add(b); padj[b].add(a)
    nid = max(X) + 1
    eid = max(int(v) for v in re.findall(r'<elem id="(\d+)"', t)) + 1
    newX, ties, new_anchor, blocks, edata, report = {}, '', [], '', '', []
    surf_add, surf_drop = defaultdict(list), set()
    surfs = {mm.group(1): [frozenset(map(int, f.split(','))) for f in re.findall(r'>([\d,]+)</(?:tri3|quad4)>', mm.group(2))]
             for mm in re.finditer(r'<Surface name="([^"]+)">(.*?)</Surface>', t, re.S)}
    ad = re.findall(r'<ElementData type="mat_axis" elem_set="([^"]+)">(.*?)</ElementData>', t, re.S)
    for s in STRAPS:
        W = {k: c for k, (n, c) in E.items() if n == s}
        topof = {c[i]: c[i + 3] for c in W.values() for i in range(3)}
        be = Counter(frozenset((c[i], c[(i + 1) % 3])) for c in W.values() for i in range(3))
        badj = defaultdict(set)
        for e, n in be.items():
            if n == 1:
                a, b = tuple(e); badj[a].add(b); badj[b].add(a)
        # the two boundary paths curve node -> anchor through free nodes; medial = the one whose anchor end is nearer x = 0
        paths = []
        for c0 in [v for v in badj if iscurve(v)]:
            for nx in badj[c0]:
                if iscurve(nx) or nx in anchors:
                    continue
                p = [c0, nx]
                while p[-1] not in anchors and not iscurve(p[-1]):
                    p.append(next(v for v in badj[p[-1]] if v != p[-2]))
                if p[-1] in anchors:          # a path back onto the curve is the curve edge itself
                    paths.append(p)
        med = min(paths, key=lambda p: abs(X[p[-1]][0]))
        C, A = med[0], med[-1]
        # curve edge: rim from C away from the strap to the first AVW node
        cp = [C]
        while '_PickedSet347' not in dn[cp[-1]]:
            nx = [v for v in rim[cp[-1]] if v not in cp and s not in dn[v]]
            assert nx, (s, cp)
            cp.append(min(nx, key=lambda v: abs(X[v][0])))
        # arch edge: along the PM_Plane boundary from the point nearest A, away from the other anchors, ARCH mm
        b0 = min(padj, key=lambda v: np.linalg.norm(X[v] - X[A]))
        prev_anc = med[-2]
        other = [v for v in anchors if s in dn[v] and v != A]
        oc = np.mean([X[v] for v in other], 0)
        nb = max(padj[b0], key=lambda v: np.linalg.norm(X[v] - oc))
        if np.linalg.norm(X[b0] - oc) < np.linalg.norm(X[A] - oc):     # b0 lies behind A (anchor side): start at nb
            ap, L, pv, cv = [X[A], X[nb]], np.linalg.norm(X[nb] - X[A]), b0, nb
        else:
            ap, L, pv, cv = [X[A], X[b0], X[nb]], np.linalg.norm(X[b0] - X[A]) + np.linalg.norm(X[nb] - X[b0]), b0, nb
        while L < arch:
            nx = next(v for v in padj[cv] if v != pv)
            L += np.linalg.norm(X[nx] - X[cv]); ap.append(X[nx]); pv, cv = cv, nx
        seg = np.linalg.norm(np.diff(np.array(ap), axis=0), axis=1); cum = np.r_[0, np.cumsum(seg)]
        i = np.searchsorted(cum, arch) - 1
        ap = ap[:i + 1] + [ap[i] + (arch - cum[i]) / seg[i] * (ap[i + 1] - ap[i])]
        ns = len(med) - 1
        Cs = resample([X[v] for v in cp], nw)
        As = resample(ap, nw)
        M = [X[v] for v in med]
        smed = np.r_[0, np.cumsum(np.linalg.norm(np.diff(M, axis=0), axis=1))]; smed /= smed[-1]
        F0, F1 = X[cp[-1]], As[-1][0]
        B = {}
        for i in range(ns + 1):
            for j in range(nw + 1):
                u, sv = j / nw, smed[i]
                if j == 0:
                    B[i, j] = med[i]; continue
                if i == 0:
                    p, si, f = Cs[j]
                    a, b = cp[si], cp[si + 1]
                    if f < 1e-6 or f > 1 - 1e-6:
                        B[i, j] = a if f < 0.5 else b; continue
                    newX[nid] = p; ties += lc(nid, [(a, 1 - f), (b, f)]); B[i, j] = nid; nid += 1; continue
                if i == ns:
                    newX[nid] = As[j][0]; new_anchor.append(nid); B[i, j] = nid; nid += 1; continue
                Fe = (1 - sv) * F0 + sv * F1
                p = ((1 - u) * M[i] + u * Fe + (1 - sv) * Cs[j][0] + sv * As[j][0]
                     - ((1 - u) * (1 - sv) * M[0] + u * (1 - sv) * F0 + (1 - u) * sv * M[-1] + u * sv * F1))
                newX[nid] = p; B[i, j] = nid; nid += 1
        Tn = {}
        P = lambda v: newX[v] if v in newX else X[v]   # noqa: E731
        for i in range(ns + 1):
            th = X[topof[med[i]]] - X[med[i]]
            for j in range(nw + 1):
                if j == 0:
                    Tn[i, j] = topof[med[i]]
                else:
                    newX[nid] = P(B[i, j]) + th; Tn[i, j] = nid; nid += 1
        # wedges
        wed = []
        for i in range(ns):
            for j in range(nw):
                q = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                for tri in ((q[0], q[1], q[2]), (q[0], q[2], q[3])):
                    c = [B[k] for k in tri] + [Tn[k] for k in tri]
                    if wvol(np.array([P(v) for v in c])) < 0:
                        c = [c[0], c[2], c[1], c[3], c[5], c[4]]
                    assert wvol(np.array([P(v) for v in c])) > 0
                    wed.append((eid, c)); eid += 1
        # fibre directions from the nearest strap element
        par = dict(ad)[s]
        rows = re.findall(r'<elem lid="(\d+)">(.*?)</elem>', par)
        wl = list(W.items())
        cen = np.array([np.mean([X[v] for v in c], 0) for _, c in wl])
        nd = s.replace('_solid', '_wide')
        edata += f'\t\t<ElementData type="mat_axis" elem_set="{nd}">\n'
        for lid, (k, c) in enumerate(wed, 1):
            near = int(np.argmin(np.linalg.norm(cen - np.mean([P(v) for v in c], 0), axis=1)))
            edata += f'\t\t\t<elem lid="{lid}">{rows[near][1]}</elem>\n'
        edata += '\t\t</ElementData>\n'
        blocks += f'\t\t<Elements type="penta6" name="{nd}">\n' + ''.join(
            f'\t\t\t<elem id="{k}">' + ','.join(map(str, c)) + '</elem>\n' for k, c in wed) + '\t\t</Elements>\n'
        # surfaces
        bt = {frozenset(c[:3]) for c in W.values()}; tt = {frozenset(c[3:]) for c in W.values()}
        sq = {frozenset((c[a], c[b], c[b + 3], c[a + 3])) for c in W.values() for a, b in ((0, 1), (1, 2), (2, 0))}
        medsq = {frozenset((med[i], med[i + 1], topof[med[i + 1]], topof[med[i]])) for i in range(ns)}
        far = [frozenset((B[i, nw], B[i + 1, nw], Tn[i + 1, nw], Tn[i, nw])) for i in range(ns)]
        for name, fs in surfs.items():
            fs = set(fs)
            hb, ht, hs = bool(fs & bt), bool(fs & tt), bool(fs & (sq - medsq))
            for k, c in wed:
                if hb:
                    surf_add[name].append(('tri3', [c[0], c[2], c[1]]))
                if ht:
                    surf_add[name].append(('tri3', c[3:]))
            if hs:
                for f in far:
                    k, c = next((k, c) for k, c in wed if f <= set(c))
                    lat = next(((c[a], c[b], c[b + 3], c[a + 3]) for a, b in ((0, 1), (1, 2), (2, 0))
                                if frozenset((c[a], c[b], c[b + 3], c[a + 3])) == f))
                    surf_add[name].append(('quad4', list(lat)))
            for f in fs & medsq:
                surf_drop.add((name, f))
        report.append(f'{s}: medial edge {len(med)} nodes ({C} -> {A}); curve edge {cp} ({np.linalg.norm(X[cp[-1]] - X[C]):.2f} mm '
                      f'chord); arch edge {arch} mm from anchor {A} (end {np.round(As[-1][0], 2)}); {len(wed)} wedges, '
                      f'vol min {min(wvol(np.array([P(v) for v in c])) for _, c in wed):.4f} mm3')
    # ---- write
    b = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:b] + '\t\t<Nodes name="strap_wide_nodes">\n' + ''.join(
        f'\t\t\t<node id="{v}">' + ','.join(f'{x:.7g}' for x in newX[v]) + '</node>\n' for v in sorted(newX)) + '\t\t</Nodes>\n' + t[b:]
    me = t.index('</Mesh>'); b = t.rindex('</Elements>', 0, me) + len('</Elements>\n')
    t = t[:b] + blocks + t[b:]
    for s in STRAPS:
        mm = re.search(r'\t\t<SolidDomain name="%s" (mat="[^"]+") />\n' % s, t)
        t = t[:mm.end()] + f'\t\t<SolidDomain name="{s.replace("_solid", "_wide")}" {mm.group(1)} />\n' + t[mm.end():]
    b = t.rindex('</ElementData>') + len('</ElementData>\n')
    t = t[:b] + edata + t[b:]
    t = re.sub(r'(<NodeSet name="curve_end_PM_anchors">)([^<]*)(</NodeSet>)',
               lambda mm: mm.group(1) + mm.group(2).rstrip() + ',' + ','.join(map(str, new_anchor)) + mm.group(3), t)

    def surf(mm):
        name = mm.group(1); out, n = [], 0
        for fm in re.finditer(r'<(quad4|tri3) id="\d+">([\d,]+)</\1>', mm.group(2)):
            if (name, frozenset(map(int, fm.group(2).split(',')))) in surf_drop:
                continue
            n += 1; out.append(f'\t\t\t<{fm.group(1)} id="{n}">{fm.group(2)}</{fm.group(1)}>\n')
        for kind, c in surf_add.get(name, []):
            n += 1; out.append(f'\t\t\t<{kind} id="{n}">' + ','.join(map(str, c)) + f'</{kind}>\n')
        return mm.group(0)[:mm.start(2) - mm.start(0)] + ''.join(out) + '\t\t</Surface>'
    t = re.sub(r'<Surface name="([^"]+)">\n(.*?)\t\t</Surface>', surf, t, flags=re.S)
    # the tie block copies strap_refine_ties' settings; build_batch443 (strap merged onto the fine rim, L440g on) removed that
    # block, so then spring_patch_ties' (2026-10-23, L451)
    for tmpl in ('strap_refine_ties', 'spring_patch_ties'):
        blk = re.search(r'\t\t<constraint name="%s"[^>]*>\n(.*?)(?=\t\t\t<linear_constraint>)' % tmpl, t, re.S)
        if blk:
            break
    end = t.index('</constraint>\n', blk.start()) + len('</constraint>\n')
    t = t[:end] + blk.group(0).replace(tmpl, 'strap_wide_ties') + ties + '\t\t</constraint>\n' + t[end:]
    msg = ('NOT IN SOURCE (the strap is ours): fan straps widened on their medial side (build_batch413, the user\'s red area): '
           + '; '.join(report) + f'; {len(newX)} new nodes, {len(new_anchor)} new fixed anchors, '
           f'{ties.count("<linear_constraint>") // 3} curve-edge nodes tied onto the rim; surfaces: '
           + ', '.join(f'{k} +{len(v)}' for k, v in surf_add.items()) + f'; {len(surf_drop)} inner side quads removed')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch413.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    ar = float(a[a.index('--arch') + 1]) if '--arch' in a else 5.0
    nw = int(a[a.index('--nw') + 1]) if '--nw' in a else 5
    main(a[0], a[1], ar, nw)
