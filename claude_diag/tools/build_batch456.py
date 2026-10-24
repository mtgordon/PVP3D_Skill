"""Batch 456 (2026-10-23; the user on L454's springs: "some of them are attached at the wrong places ... at least one of them crosses
over another connector"; "it's the [AVW-Para] ones that attached where the mesh size was refined"; the strap springs' curve ends ->
"Outer face of the curve"). NOT IN SOURCE (spring ends):
1. PM_avw_strap_left / right_conn (build_batch454): each curve end moves to the nearest node on the curve's OUTER surface (a curve
   boundary node, the one farthest from the local canal centre among those within RAD mm), then the springs are re-paired in order
   along the edge and along the PM-arch anchors (both sorted by their position along the edge's main direction), so none cross.
2. AVW-Para-L / R_conn ends that sat on spring patches in REF (build_batch288; 14 a side, rest z -39..-15; build_batch452 then
   moved the 3 lowest onto the AVW top): the patches' tie rings were made on the coarse curve, which was refined after (build_batch412);
   each patch anchor sits exactly on the curve node of the attachment line (one element up from the curve's middle, toward the AVW;
   the user), so the spring goes straight onto that node, no patch. The single PM-side hub node stays (as in the source).
usage: py -3.10 build_batch456.py SRC_RUN NEW_RUN [--ref L451_cog_wog_s033_pm27_swide_sm] [--rad 3] [--para-only]
(--para-only: part 2 only, for a model that still has its strap, e.g. L461 = L451)"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]



def main(src, new, ref, rad, para_only=False):
    t, X, E = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    cur = {k: c for k, (n, c) in E.items() if n.startswith('canal_wrap') and len(c) == 8}
    cnt = Counter(frozenset(c[i] for i in f) for c in cur.values() for f in HF)
    bnd = sorted({v for c in cur.values() for f in HF if cnt[frozenset(c[i] for i in f)] == 1 for v in (c[i] for i in f)})
    canal = np.array([X[v] for n, c in E.values() if n.startswith('canal_wrap') or n in ('_PickedSet347', '_PickedSet64') for v in c])
    B = np.array([X[v] for v in bnd])

    def outer(p, taken=()):  # an unused curve boundary node within rad of p, far from the local canal centre and near p
        cen = canal[abs(canal[:, 2] - p[2]) < 1.0].mean(0)
        d = np.linalg.norm(B - p, axis=1); cand = [i for i in np.where(d < rad)[0] if bnd[i] not in taken]
        r = np.linalg.norm((B[cand] - cen)[:, :2], axis=1)
        return bnd[cand[int(np.argmax(r - 0.3 * d[cand]))]]
    report = []
    # 1. strap springs (none in a strap model: --para-only, L461 on L451)
    for side in (() if para_only else ('left', 'right')):
        m = re.search(r'(<DiscreteSet name="PM_avw_strap_%s_conn">)(.*?)(</DiscreteSet>)' % side, t, re.S)
        pr = [(int(a), int(b)) for a, b in re.findall(r'<delem>(\d+),(\d+)</delem>', m.group(2))]
        ends, taken = [], set()
        for e, a in pr:
            ends.append(outer(X[e], taken)); taken.add(ends[-1])
        anc = sorted({a for e, a in pr})
        P = np.array([X[v] for v in ends]); ax = np.linalg.svd(P - P.mean(0))[2][0]
        Q = np.array([X[v] for v in anc]); aq = ax if (Q[-1] - Q[0]) @ ax else ax
        ends = [v for _, v in sorted(zip([X[v] @ ax for v in ends], ends))]
        anc = [v for _, v in sorted(zip([X[v] @ aq for v in anc], anc))]
        pick = [anc[i] for i in np.round(np.linspace(0, len(anc) - 1, len(ends))).astype(int)]
        # orient: the pairing with the shorter total length (same order or reversed) does not cross
        if sum(np.linalg.norm(X[e] - X[a]) for e, a in zip(ends, pick)) > sum(np.linalg.norm(X[e] - X[a]) for e, a in zip(ends, pick[::-1])):
            pick = pick[::-1]
        body = '\n' + ''.join(f'\t\t\t<delem>{e},{a}</delem>\n' for e, a in zip(ends, pick)) + '\t\t'
        t = t[:m.start(2)] + body + t[m.end(2):]
        report.append(f'strap {side}: {len(ends)} ends on the curve outer face {ends}, re-paired in order')
    # 2. AVW-Para ends on spring patches (in REF, before build_batch452 moved 3 of them): each patch anchor was put at the curve
    # node of the attachment line (one element up from the curve's middle, toward the AVW; the user, 10-23), and that node still
    # exists in the refined curve -> the spring goes straight onto it (no patch)
    _, XR, ER = load(ref)
    XR = {n: np.asarray(p, float) for n, p in dict(XR).items()}
    inel_r = {v for n, c in ER.values() for v in c}
    cn = sorted({v for n, c in E.values() if n.startswith('canal_wrap') for v in c}); C = np.array([X[v] for v in cn])
    for side in ('L', 'R'):
        fam = f'AVW-Para-{side}_conn'
        m = re.search(r'(<DiscreteSet name="%s">)(.*?)(</DiscreteSet>)' % re.escape(fam), t, re.S)
        rows = re.findall(r'(\s*<delem>)(\d+),(\d+)(</delem>)', m.group(2))
        pr = [(int(a), int(b)) for _, a, b, _ in rows]
        hub = Counter([a for a, b in pr] + [b for a, b in pr]).most_common(1)[0][0]
        mr = re.search(r'<DiscreteSet name="%s">(.*?)</DiscreteSet>' % re.escape(fam), open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', ref, ref + '.feb'), encoding='utf-8').read(), re.S)
        prr = [(int(a), int(b)) for a, b in re.findall(r'<delem>(\d+),(\d+)</delem>', mr.group(1))]
        assert len(prr) == len(pr), fam
        near = [b if a == hub else a for a, b in pr]
        moved = []
        for k, (a, b) in enumerate(prr):
            rn = b if a == hub else a
            if rn in inel_r:
                continue                                                    # a curve node in REF: unchanged
            d = np.linalg.norm(C - XR[rn], axis=1); i = int(np.argmin(d))
            assert d[i] < 0.05, (fam, rn, d[i])
            moved.append((near[k], cn[i])); near[k] = cn[i]
        out = ''.join(f'{o}{hub if a == hub else near[k]},{hub if b == hub else near[k]}{c}' for k, (o, a, b, c) in
                      enumerate((o, int(a), int(b), c) for o, a, b, c in rows)) + re.search(r'\s*$', m.group(2)).group(0)
        t = t[:m.start(2)] + out + t[m.end(2):]
        report.append(f'{fam}: {len(moved)} ends from spring patches (REF {ref}) onto the attachment-line node at the patch anchor '
                      f'(old -> node) {moved}')
    msg = 'NOT IN SOURCE (spring ends): ' + '; '.join(report)
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch456.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    opt = lambda n, d: (a[a.index(n) + 1] if n in a else d)  # noqa: E731
    main(a[0], a[1], opt('--ref', 'L451_cog_wog_s033_pm27_swide_sm'), float(opt('--rad', 3)), '--para-only' in a)
