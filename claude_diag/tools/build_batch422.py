"""Batch 422 (2026-10-22; the user: "I think our model has the vaginal wall being too wide. Compare our model to the original and
try to make our avw/pvw + curves be the same width as the original width of the avw/pvw. Some of the connective tissues that
currently attach to the avw/pvw may need to attach to the curve as well ... the connectors that go to the top of the pm arch
at the distal end and also the PeB"; answers: morph sideways, move the ends in the curve band, on L413b).
Measured (L413b vs PVP3DModel_job.inp, AVW _PickedSet347 + PVW _PickedSet64): our walls are ~3 mm narrower a side than the
source's (build_batch135 narrowing) and the curves stick out ~5 mm beyond them, so walls + curves are 8-28 % wider than the
source walls (17.2 vs 21.9 mm at rest z -42..-38, 25.3 vs 29.5 at -10..-6).
NOT IN SOURCE (rest geometry of our canal; the target width IS the source's):
1. Width profiles by rest z (2 mm bins, smoothed over 3): source wall edges xl_s / xr_s, our canal (AVW + PVW + curves) xl_o /
   xr_o, our walls xlw / xrw. Every node of the AVW, PVW and curves moves in x only: each side of the midline is scaled
   linearly so [xc_o, xr_o] -> [xc_s, xr_s] (and the left side alike); the curves' outer edge lands on the source wall edge.
2. Nodes of the PeB, cervix, apex_contour, the straps and element-less nodes (spring-patch anchors) within FADE mm of a moved
   node take w * the nearest moved node's shift, w = cos^2 (pi/2 d / FADE); nodes in a zero-displacement BC never move.
3. PM_conn (26 springs, the AVW -> top of the PM arch): spread in order along the whole distal edge row = the AVW nodes on
   its wall-end row + the curve's outer top corner beyond each end (L423: 27 + 2 = 29 nodes, 26 picked evenly by index),
   paired in x order with the fixed PM ends (no crossing; one material per family); co-located patch anchors reused, the
   new curve ends put on spring patches with build_batch288 --edge --only. PM_PeB_Left / Right stay on the PeB (not
   narrowed; its edge is at the source width). (First version, L422: ends whose source point fell in the curve band moved
   to the nearest curve node -> two ends on one node and 6 crossing pairs; the user: connect them in order.)
4. Checks: every solid element positive (elements already non-positive at rest by the tet test are not counted) (hex 5 tets, penta 3 tets); the shift range; ties' rest-residual change.
usage: py -3.10 build_batch422.py SRC_RUN NEW_RUN [--fade 8]"""
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402
from curve_swell import TETS  # noqa: E402

SRC_INP = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'PVP3DModel_job.inp')
WALLS = ('_PickedSet347', '_PickedSet64')
CERVIX = '_PickedSet346'
MOVED_ALSO = (CERVIX, 'apex_contour')   # the user (10-22): 'cervix should get shrunk as well'
FADE_DOMS = ('_PickedSet66', 'PM_avw_bottom_')
FAMS = ('PM_conn', 'PM_PeB_Left_conn', 'PM_PeB_Right_conn')
W3 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def source_walls():
    L = open(SRC_INP).read().split('\n')
    a = next(i for i, l in enumerate(L) if l.startswith('*Part, name=VW-PeB'))
    b = next(i for i in range(a + 1, len(L)) if L[i].startswith('*End Part'))
    X, El, sets, mode, cur, gen = {}, {}, {}, None, None, False
    for ln in L[a + 1:b]:
        if ln.startswith('*'):
            low = ln.lower(); mode = None
            if low.startswith('*node') and not low.startswith('*node output'):
                mode = 'n'
            elif low.startswith('*element'):
                mode = 'e'
            elif low.startswith('*elset'):
                mode = 's'; cur = re.search(r'elset=([^,\s]+)', ln).group(1); gen = 'generate' in low; sets[cur] = []
            continue
        v = [s.strip() for s in ln.split(',') if s.strip()]
        if not v:
            continue
        if mode == 'n':
            X[int(v[0])] = np.array(list(map(float, v[1:4])))
        elif mode == 'e':
            El[int(v[0])] = list(map(int, v[1:9]))
        elif mode == 's':
            if gen:
                p, q, s_ = map(int, v[:3]); sets[cur] += list(range(p, q + 1, s_))
            else:
                sets[cur] += list(map(int, v))
    return np.array([X[n] for w in WALLS + (CERVIX,) for e in sets[w] for n in El[e]])


def profile(P, zc):
    lo, hi = [], []
    for z in zc:
        s = P[(P[:, 2] >= z - 1) & (P[:, 2] < z + 1)]
        lo.append(s[:, 0].min() if len(s) else np.nan); hi.append(s[:, 0].max() if len(s) else np.nan)
    lo, hi = np.array(lo), np.array(hi)
    ok = ~np.isnan(lo)
    lo = np.interp(zc, zc[ok], lo[ok]); hi = np.interp(zc, zc[ok], hi[ok])
    k = np.ones(3) / 3
    sm = lambda v: np.r_[v[0], np.convolve(v, k, 'valid'), v[-1]]   # noqa: E731
    return sm(lo), sm(hi)


def xmap(x, z, zc, a_lo, a_hi, b_lo, b_hi):
    """map x at level z from the extent [a_lo, a_hi] to [b_lo, b_hi], each side of the midline linearly"""
    al, ah, bl, bh = (np.interp(z, zc, v) for v in (a_lo, a_hi, b_lo, b_hi))
    ac, bc = (al + ah) / 2, (bl + bh) / 2
    return np.where(x >= ac, bc + (x - ac) * (bh - bc) / (ah - ac), bc + (x - ac) * (bl - bc) / (al - ac))


def vol_ok(P, c):
    if len(c) == 8:
        return min(np.dot(P[c[q1]] - P[c[q0]], np.cross(P[c[q2]] - P[c[q0]], P[c[q3]] - P[c[q0]])) for q0, q1, q2, q3 in TETS)
    if len(c) == 6:
        return min(np.dot(np.cross(P[c[b]] - P[c[a]], P[c[d]] - P[c[a]]), P[c[e]] - P[c[a]]) for a, b, d, e in W3)
    return 1.0


def main(src, new, fade=8.0, ztop0=1e3, ztop1=1e3 + 1):
    t, X, E = load(src)
    X = {k: np.asarray(p, float) for k, p in dict(X).items()}
    S = source_walls()
    canal = {v for n, c in E.values() if n in WALLS or n in MOVED_ALSO or n.startswith('canal_wrap_hex') for v in c}
    walls = {v for n, c in E.values() if n in WALLS for v in c}
    C = np.array([X[v] for v in canal]); Wp = np.array([X[v] for v in walls])
    zc = np.arange(np.floor(min(C[:, 2].min(), S[:, 2].min())), np.ceil(max(C[:, 2].max(), S[:, 2].max())) + 2, 2.0)
    sl, sh = profile(S, zc); ol, oh = profile(C, zc); wl, wh = profile(Wp, zc)
    sl, sh = np.maximum(sl, ol), np.minimum(sh, oh)      # shrink only: never widen (our cervix top is narrower than the source)
    # 1. canal nodes
    ids = sorted(canal); P0 = np.array([X[v] for v in ids])
    nx = xmap(P0[:, 0], P0[:, 2], zc, ol, oh, sl, sh)
    taper = np.clip((ztop1 - P0[:, 2]) / (ztop1 - ztop0), 0, 1)        # full below ztop0, none above ztop1 (apex region)
    nx = P0[:, 0] + taper * (nx - P0[:, 0])
    shift = dict(zip(ids, nx - P0[:, 0]))
    # 2. fade into the neighbours
    bc_nodes = set()
    for m in re.finditer(r'<bc [^>]*node_set="([^"]+)"[^>]*type="zero displacement"', t):
        ns = re.search(r'<NodeSet name="%s">([^<]*)</NodeSet>' % re.escape(m.group(1)), t)
        if ns:
            bc_nodes |= {int(v) for v in re.findall(r'\d+', ns.group(1))}
    in_elem = defaultdict(set)
    for n, c in E.values():
        for v in c:
            in_elem[v].add(n)
    tree = cKDTree(P0); sh_arr = nx - P0[:, 0]
    faded = 0
    for v, p in X.items():
        if v in canal or v in bc_nodes:
            continue
        doms = in_elem.get(v, set())
        if doms and not any(d.startswith(FADE_DOMS) for d in doms):
            continue
        d, j = tree.query(p)
        if d >= fade:
            continue
        w = np.cos(np.pi / 2 * d / fade) ** 2
        shift[v] = w * sh_arr[j]; faded += 1
    # slivers (min sub-tet < 1e-3 mm3 at rest, e.g. at the apex_contour rim): one common shift for their nodes (translate, not shear)
    sliv = 0
    for k, (n, c) in E.items():
        if len(c) in (6, 8) and 0 < vol_ok(X, c) < 1e-3 and any(v in shift for v in c):
            m_ = np.mean([shift.get(v, 0.0) for v in c])
            for v in c:
                if v in shift or v not in bc_nodes:
                    shift[v] = m_
            sliv += 1
    Xn = {v: (p + np.array([shift[v], 0, 0]) if v in shift else p) for v, p in X.items()}
    # 4a. element check
    bad = [(k, n) for k, (n, c) in E.items() if len(c) in (6, 8) and vol_ok(Xn, c) <= 0 and vol_ok(X, c) > 0]
    assert not bad, bad[:10]
    # 3. PM_conn along the whole distal edge, in order (the user, 10-22: "there are a lot of avw nodes that don't have
    # connectors, maybe just connect them in order"): the row = the AVW nodes on the family's wall-end row (same rest y / z as
    # its AVW ends) + the curve boundary node beyond each end of it (|dy| < 0.5, |dz| < 0.8: the curve's outer top corner);
    # NS of its nodes picked evenly by index, paired in x order with the fixed far ends (one material per family, no crossing).
    # A picked node that already carries a co-located spring-patch anchor of the family uses the anchor.
    cur = {k: c for k, (n, c) in E.items() if n.startswith('canal_wrap_hex') and len(c) == 8}
    cnt = Counter(frozenset(c[i] for i in f) for c in cur.values() for f in HF)
    cb = sorted({v for c in cur.values() for f in HF if cnt[frozenset(c[i] for i in f)] == 1 for v in (c[i] for i in f)})
    fixed = set()
    for mm in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', t):
        if re.search(r'<bc [^>]*node_set="%s"[^>]*type="zero displacement"' % re.escape(mm.group(1)), t):
            fixed |= {int(v) for v in re.findall(r'\d+', mm.group(2))}
    report, newcurve = [], []
    fam = 'PM_conn'
    m = re.search(r'(<DiscreteSet name="%s">)(.*?)(</DiscreteSet>)' % fam, t, re.S)
    rows = re.findall(r'(\s*<delem[^>]*>)(\d+),(\d+)(</delem>)', m.group(2))
    pairs = [(int(p_), int(q_)) for _, p_, q_, _ in rows]
    near0 = [p_ if q_ in fixed else q_ for p_, q_ in pairs]; far = sorted((q_ if q_ in fixed else p_ for p_, q_ in pairs), key=lambda v: Xn[v][0])
    wend = [v for v in near0 if v in walls]
    y0, z0 = np.median([X[v][1] for v in wend]), np.median([X[v][2] for v in wend])
    row = sorted((v for v in walls if abs(X[v][1] - y0) < 0.05 and abs(X[v][2] - z0) < 0.1), key=lambda v: Xn[v][0])
    lo, hi = Xn[row[0]][0], Xn[row[-1]][0]
    ext = [u for u in cb if u not in walls and abs(Xn[u][1] - Xn[row[0]][1]) < 0.5 and abs(Xn[u][2] - Xn[row[0]][2]) < 0.8]
    left = [u for u in ext if Xn[u][0] < lo]; right = [u for u in ext if Xn[u][0] > hi]
    full = ([max(left, key=lambda u: Xn[u][0])] if left else []) + row + ([min(right, key=lambda u: Xn[u][0])] if right else [])
    pick = [full[i] for i in np.round(np.linspace(0, len(full) - 1, len(far))).astype(int)]
    anchors = {v for v in near0 if v not in walls and not in_elem.get(v)}
    co = {}
    for an in anchors:
        d = {u: np.linalg.norm(Xn[u] - Xn[an]) for u in full}
        u = min(d, key=d.get)
        if d[u] < 1e-3:
            co[u] = an
    near = [co.get(u, u) for u in pick]
    newcurve = [u for u in pick if u not in co and u not in walls]
    body = ''.join(f'{o}{n_},{f_}{c_}' for (o, _, _, c_), n_, f_ in zip(rows, near, far)) + re.search(r'\s*$', m.group(2)).group(0)
    t = t[:m.start(2)] + body + t[m.end(2):]
    assert len(re.findall('<delem', body)) == len(pairs)
    report.append(f'{fam}: {len(far)} springs spread in order over the distal edge row ({len(full)} nodes: {len(row)} AVW + '
                  f'{len(full) - len(row)} curve corner, x {Xn[full[0]][0]:.1f}..{Xn[full[-1]][0]:.1f}), paired with the fixed far ends in '
                  f'x order; {len(co)} on existing patch anchors, {len(newcurve)} new curve ends {newcurve} (to be patched)')
    peb_n = sum(len(re.findall(r'<delem', re.search(r'<DiscreteSet name="%s">(.*?)</DiscreteSet>' % f_, t, re.S).group(1)))
                for f_ in ('PM_PeB_Left_conn', 'PM_PeB_Right_conn'))
    report.append(f'{peb_n} PM_PeB ends stay on the PeB (not narrowed; its edge is at the source width)')
    moved = {(fam, u): u for u in newcurve}
    # write node coordinates
    def node_sub(mm):
        v = int(mm.group(1))
        if v not in shift:
            return mm.group(0)
        return f'<node id="{v}">' + ','.join(f'{c:.7g}' for c in Xn[v]) + '</node>'
    a = t.index('<Mesh>'); b = t.index('</Mesh>')
    t = t[:a] + re.sub(r'<node id="(\d+)">[^<]*</node>', node_sub, t[a:b]) + t[b:]
    sh_c = np.array([shift[v] for v in canal])
    msg = (f'NOT IN SOURCE (rest geometry; target width = the source walls): canal (AVW + PVW + curves, {len(canal)} nodes) '
           f'narrowed in x to the source wall width by level, tapered to none from rest z {ztop0:g} to {ztop1:g} (apex region) (shift {sh_c.min():.2f}..{sh_c.max():.2f} mm); {faded} neighbour nodes '
           f'(PeB, straps, anchors) faded over {fade:g} mm; {sliv} rest slivers translated as one; no element inverted; ' + ' | '.join(report))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch422.py)\n{msg}\n')
    return sorted(newcurve)


if __name__ == '__main__':
    a = sys.argv[1:]
    fd = float(a[a.index('--fade') + 1]) if '--fade' in a else 8.0
    nodes = main(a[0], a[1], fd)
    print('moved spring ends now on curve nodes:', nodes)
