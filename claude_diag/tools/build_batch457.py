"""Batch 457 (2026-10-23; the user: "some of the connectors from the top still attach at nodes that are not on the distal face. Just
go to each node on the outside of the distal face and put a connector that goes to the arch"). NOT IN SOURCE (spring layout):
PM_conn (AVW top row -> PM arch, build_batch448) and the strap springs PM_avw_strap_left / right_conn (build_batch454 / 456) are
replaced by ONE family, PM_conn: a spring from every node on the outer edge of the canal's distal face (the AVW's outer distal row,
its two seam nodes with the curves included, + each curve's outer distal rim up to the PVW, i.e. the curve-only nodes of the distal
end faces that lie on the outward side surface; the PVW's own end, on the PeB, left out) to the PM arch. The arch points are the old
PM_conn far ends + the strap springs' arch anchors; nodes and arch points are both sorted by angle around the canal's distal centre
and paired in order (arch points picked evenly by index), so no two cross. Material PM_conn_mat (table unchanged), <scale> = the
old family totals (PM_conn: N x scale; strap springs: N x scale, both sides) / the new count, so the total pull is unchanged.
usage: py -3.10 build_batch457.py SRC_RUN NEW_RUN [--zend -33]"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
AVW, PVW = '_PickedSet347', '_PickedSet64'


def main(src, new, zend):
    t, X, E = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    dom = {}
    for n, c in E.values():
        for v in c:
            dom.setdefault(v, set()).add(n)
    els = [(n, c) for n, c in E.values() if (n.startswith('canal_wrap') or n in (AVW, PVW)) and len(c) == 8]
    cnt = Counter(frozenset(c[i] for i in f) for n, c in els for f in HF)
    faces = []
    for n, c in els:
        cen = np.mean([X[v] for v in c], 0)
        for f in HF:
            if cnt[frozenset(c[i] for i in f)] == 1:
                P = np.array([X[c[i]] for i in f]); nr = np.cross(P[2] - P[0], P[3] - P[1]); nr /= np.linalg.norm(nr)
                if (P.mean(0) - cen) @ nr < 0:
                    nr = -nr
                faces.append((tuple(c[i] for i in f), P.mean(0), nr))
    end = [F for F in faces if F[1][2] < zend and F[2][2] < -0.7]                      # distal end faces
    endn = {v for F in end for v in F[0]}
    cen = np.array([X[v] for v in endn]).mean(0)
    outer = set()
    for F in faces:
        if F not in end and F[1][2] < zend + 3 and (F[1][:2] - cen[:2]) @ F[2][:2] > 0:   # outward side faces
            outer |= {v for v in F[0] if v in endn}
    avw_row_y = max(X[v][1] for v in endn if dom[v] == {AVW})
    keep = [v for v in outer if (AVW in dom[v] and abs(X[v][1] - avw_row_y) < 0.06)        # AVW outer row (+ seam nodes)
            or (AVW not in dom[v] and PVW not in dom[v])]                                 # curve-only rim nodes
    ang = lambda p: np.arctan2(p[1] - cen[1], p[0] - cen[0])  # noqa: E731
    # angles measured from the PVW side (the gap) so the loop opens there
    gap = ang(np.array([X[v] for v in endn if dom[v] == {PVW}]).mean(0))
    rel = lambda p: (ang(p) - gap) % (2 * np.pi)  # noqa: E731
    keep.sort(key=lambda v: rel(X[v]))
    # arch points: PM_conn far ends + strap spring anchors
    fixed = set()   # fully fixed only (x, y and z): BC-VW-mid pins the wall's midline nodes in some directions only
    for mm in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', t):
        bc = re.search(r'<bc [^>]*node_set="%s"[^>]*type="zero displacement">(.*?)</bc>' % re.escape(mm.group(1)), t, re.S)
        if bc and all(re.search(r'<%s_dof>1</%s_dof>' % (d, d), bc.group(1)) for d in 'xyz'):
            fixed |= {int(v) for v in re.findall(r'\d+', mm.group(2))}
    tot, arch, old = 0.0, set(), []
    for fam in ('PM_conn', 'PM_avw_strap_left_conn', 'PM_avw_strap_right_conn'):
        m = re.search(r'<DiscreteSet name="%s">(.*?)</DiscreteSet>' % fam, t, re.S)
        pr = [(int(a), int(b)) for a, b in re.findall(r'<delem>(\d+),(\d+)</delem>', m.group(1))]
        arch |= {v for p in pr for v in p if v in fixed}
        dm = re.search(r'<discrete dmat="(\d+)" discrete_set="%s" />' % fam, t).group(1)
        mat = re.findall(r'<discrete_material id="\d+"[^>]*>.*?</discrete_material>', t, re.S)[int(dm) - 1]
        sc = float(re.search(r'<scale>([^<]+)</scale>', mat).group(1)); tot += len(pr) * sc
        old.append(f'{fam} {len(pr)} x {sc:g}')
    arch = sorted(arch, key=lambda v: rel(X[v]))
    pick = [arch[i] for i in np.round(np.linspace(0, len(arch) - 1, len(keep))).astype(int)]

    def gap_(p1, p2, q1, q2):  # closest approach of two segments (sampled)
        return min(np.linalg.norm(p1 + s * (p2 - p1) - (q1 + np.clip(((p1 + s * (p2 - p1) - q1) @ (q2 - q1)) / ((q2 - q1) @ (q2 - q1)), 0, 1)
                                                         * (q2 - q1))) for s in np.linspace(0, 1, 81))
    swaps = []
    for i in range(len(keep) - 1):  # neighbours whose arch points sit out of order (near-equal angles): swap if they nearly touch
        if pick[i] != pick[i + 1] and gap_(X[keep[i]], X[pick[i]], X[keep[i + 1]], X[pick[i + 1]]) < 0.1:
            pick[i], pick[i + 1] = pick[i + 1], pick[i]; swaps.append(keep[i])
    scale = tot / len(keep)
    # write: PM_conn gets the new springs and scale; the strap spring sets and their discrete lines go
    body = '\n' + ''.join(f'\t\t\t<delem>{e},{a}</delem>\n' for e, a in zip(keep, pick)) + '\t\t'
    m = re.search(r'(<DiscreteSet name="PM_conn">)(.*?)(</DiscreteSet>)', t, re.S)
    t = t[:m.start(2)] + body + t[m.end(2):]
    pmc = re.search(r'<discrete_material id="\d+" name="PM_conn_mat"[^>]*>.*?</discrete_material>', t, re.S)
    t = t[:pmc.start()] + re.sub(r'<scale>[^<]+</scale>', f'<scale>{scale:.6g}</scale>', pmc.group(0)) + t[pmc.end():]
    for fam in ('PM_avw_strap_left_conn', 'PM_avw_strap_right_conn'):
        t, n1 = re.subn(r'\t\t<DiscreteSet name="%s">.*?</DiscreteSet>\n' % fam, '', t, flags=re.S)
        t, n2 = re.subn(r'\t\t<discrete dmat="\d+" discrete_set="%s" />\n' % fam, '', t)
        assert (n1, n2) == (1, 1), fam
    na = sum(1 for v in keep if AVW in dom[v])
    msg = (f'NOT IN SOURCE (spring layout): the top connectors as one family PM_conn: {len(keep)} springs, one from every node on the '
           f'outer edge of the canal\'s distal face ({na} on the AVW outer distal row incl. the 2 seam nodes, {len(keep) - na} on the curves\' '
           f'outer distal rims; the PVW end left out) to the PM arch ({len(arch)} arch points: the old PM_conn far ends + strap-spring '
           f'anchors), paired in order around the canal (no crossing); replaces {", ".join(old)}; PM_conn_mat scale {scale:.4g} '
           f'(total pull unchanged); the strap spring sets removed; {len(swaps)} neighbour pair(s) swapped at the arch to avoid a crossing. One change.')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch457.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--zend') + 1]) if '--zend' in a else -33)
