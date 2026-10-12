"""Batch 166 (2026-10-05, the user, with a screenshot of the right curve's distal end: "the AVW curve doesn't actually reach
the end. Please redo it to make it reach the end. Also, add attachments so that every node on the curve (except for the ones
that attach to the AVW and PVW) have attachments"; then: change the PM_avw_bottom springs' anchors "to be the nearest point on
the outside of the pm arch and do the same for the new ones", the new springs "the same material properties as the
pm_conn_mat", only the outer-surface nodes "on the distal end that don't already have connectors" (~5 a side), on the
damped x0.35-wall version).
  1. geometry: build_batch143 rerun as L166_thin275_wrapend_full with the right side's distal end seam node 1720 added to
     the seam list. build_batch90.seam_nodes keeps only the two corners of each end row; at the perineal-body end that kept
     1773 on the left but dropped its mirror 1720 on the right, so the right curve (55 stations against 56) stopped one
     station short of the AVW's end (1720 / 1771 / the AVW top node had no curve). Everything else as L143.
  2. PM_conn on the PM's outer arc, as L162 (build_batch149.outer_arc on step 1's stage file before its PM_conn step, then
     step 1's lumen surfaces, which include the curves' inner faces).
  3. the 4 PM_avw_bottom springs (source CONN3D2 246-249, the PM-LA law at full, as L162_..._full): curve ends as L163 (the
     upper spring of each side to the outer-surface node on the end face of the curve's AVW half nearest its PM end, the
     lower to the PVW half's outer-surface node nearest its PM end); their PM ends moved to the nearest point on the PM's
     outer arc (the clamped edge) to their curve end and fixed there; their PM_anchor_ties entries dropped.
  4. new springs (set curve_end_PM_conn, material PM_conn_mat, as PM_conn): from every outer-surface node of the distal end
     face of each curve that is not on the AVW / PVW (in a wall, or tied / embedded to one) and has no spring, to a new node
     at the nearest point on the PM's outer arc, fixed.
  5. self-contact of the curves' outer surface (build_batch164.add_self_contact, recomputed on this geometry).
  6. wall damping C 200/s on the walls (build_batch165.wall_damping), walls x0.35 (build_batch145), healthy / P1 / P2.
  L166_thin275_wrapend_full                       step 1 (as L143, PM_conn inner arc)
  L166_wrapfull_outer                             + step 2
  L166_wrapfull_outer_endspr                      + steps 3-4
  L166_wrapfull_endspr_self_damp                  + step 5 + the wall damping (source wall stiffness)
  L166_wrapfull_endspr_walls035_{healthy,P1,P2}   + walls x0.35 (+ the cases)
usage: py -3.10 build_batch166.py"""
import os
import re
import shutil
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
import build_batch90  # noqa: E402
import build_batch135 as b135  # noqa: E402
import build_batch136 as b136  # noqa: E402
from build_batch113 import pm_arcs  # noqa: E402
from build_batch149 import outer_arc  # noqa: E402
from build_batch162 import write, copy_surface, elems, block  # noqa: E402
from build_batch164 import add_self_contact, outer_faces, HF, CURVES, WALLS, constrained  # noqa: E402
from build_batch165 import wall_damping  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
TOOLS = os.path.dirname(os.path.abspath(__file__))
GEOM = 'L166_thin275_wrapend_full'
OUTER = 'L166_wrapfull_outer'
ENDSPR = 'L166_wrapfull_outer_endspr'
SELF = 'L166_wrapfull_endspr_self_damp'
ADD_SEAM = (1720,)
PMAB = ('PM_avw_bottom_left_conn', 'PM_avw_bottom_right_conn')
R_INNER = 2.67            # the curves' inner radius (build_batch135 log, median)
_seam_nodes = build_batch90.seam_nodes


def seam_fixed(m, **kw):
    seam, fP, X = _seam_nodes(m, **kw)
    return list(seam) + [n for n in ADD_SEAM if n not in seam], fP, X


def step1():
    for nm in (GEOM, GEOM + '_tmppre', GEOM + '_tmpmid', GEOM + '_tmpavg'):
        assert not os.path.exists(os.path.join(RUNS, nm)) and not os.path.exists(os.path.join(SCRATCH, nm)), nm
    b135.seam_nodes = seam_fixed
    b136.run_main(os.path.join(TOOLS, 'build_batch143.py'),
                  {'NAME': GEOM, 'PRE': GEOM + '_tmppre', 'MID': GEOM + '_tmpmid', 'AVG': GEOM + '_tmpavg',
                   'seam_nodes': seam_fixed})
    with open(os.path.join(RUNS, GEOM, GEOM + '.feb.changes.txt'), 'a', encoding='utf-8') as fo:
        fo.write(f'(tools/build_batch166.py) as L143_thin275_narrow3_wrapend_noflare but with node(s) {ADD_SEAM} added to the '
                 f'seam list (the right side\'s distal end seam node, dropped by seam_nodes\' end-row rule) so the right curve '
                 f'reaches the AVW\'s end\n')


def step2():
    pre = os.path.join(SCRATCH, GEOM + '_tmppre', GEOM + '_tmppre.feb')
    t0 = open(pre, encoding='utf-8').read()
    tg = open(os.path.join(RUNS, GEOM, GEOM + '.feb'), encoding='utf-8').read()
    t, msg = outer_arc(t0, pre)
    for s in ('SlidingElastic1Primary', 'SlidingElastic1Secondary'):
        t = copy_surface(t, tg, s)
    write(OUTER, t, f'{OUTER} (tools/build_batch166.py) from {os.path.relpath(pre, JOBS)} (= {GEOM} before its PM_conn step):\n'
                    f'{open(pre + ".changes.txt", encoding="utf-8").read()}{msg}\nthe curves\' inner faces in the lumen '
                    f'contact: SlidingElastic1Primary / Secondary as in {GEOM}\n')


def arc_point(O, p):
    """nearest point to p on the polyline O"""
    best = None
    for k in range(len(O) - 1):
        a, b = O[k], O[k + 1]
        u = float(np.clip((p - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1))
        q = a + u * (b - a)
        d = np.linalg.norm(q - p)
        if best is None or d < best[0]:
            best = (d, q)
    return best[1]


def step34():
    p = os.path.join(RUNS, OUTER, OUTER + '.feb')
    t = open(p, encoding='utf-8').read()
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    _, outer, _, _ = pm_arcs(p)
    O = np.array([X[n] for n in outer])
    # the curves' surfaces
    E = {d: elems(t, d) for d in CURVES + WALLS}
    faces, _ = outer_faces(t)
    onodes = {n for f in faces for n in f}
    mem = {d: {n for e in E[d] for n in e} for d in CURVES}
    wall = {n for d in WALLS for e in E[d] for n in e}
    tied = constrained(t, 'canal_wrap_ends') | constrained(t, 'canal_wrap_embed')
    # the outer ring of the curves' distal ends: a node on the curve's end face and on its outer or lumen surface lies in at
    # most 2 curve hexes (outer-surface nodes along the curve lie in 4, through-thickness ones in 4 or more); of those at the
    # distal end (z < -30) that are not on a wall or tied / embedded to one, the outer ones are farther than the mid-radius
    # from the ring centre (lumen-side ones lie at the inner radius, ~2.7 mm). The centre: medially from the inner apex (the
    # end node in both halves nearest the midline) by the inner radius, along the apex column.
    hexcount = Counter(n for d in CURVES for e in E[d] for n in set(e))
    outer_end, rad = set(), {}
    for sgn in (-1, 1):
        ring = [n for n, k in hexcount.items() if k <= 2 and X[n][2] < -30 and np.sign(X[n][0] - 1.4) == sgn
                and n not in wall and n not in tied]
        col = sorted((n for n in ring if n in mem[CURVES[0]] and n in mem[CURVES[1]]), key=lambda n: sgn * X[n][0])
        inner_apex, outer_apex = X[col[0]], X[col[-1]]
        d = (outer_apex - inner_apex) / np.linalg.norm(outer_apex - inner_apex)
        c = inner_apex - d * R_INNER
        R = np.linalg.norm(outer_apex - c)
        for n in ring:
            rad[n] = (np.linalg.norm(X[n] - c), R)
            if rad[n][0] > 0.5 * (R_INNER + R):
                outer_end.add(n)
    # 3. PM_avw_bottom: curve ends, then anchors to the outer arc
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    blk = t[i0:i1]
    n0 = blk.count('<linear_constraint>')
    rep3, anchors, targets = [], {}, set()
    for nm in PMAB:
        i, j = block(t, 'DiscreteSet', nm)
        sb = t[i:j]
        pairs = [tuple(map(int, d.split(','))) for d in re.findall(r'<delem>([^<]+)</delem>', sb)]
        assert len(pairs) == 2
        for (anc, old), (half, pool) in zip(pairs, (('AVW', outer_end - mem[CURVES[1]]),
                                                    ('PVW', outer_end - mem[CURVES[0]]))):
            pool = {n for n in pool if n not in wall and n not in tied and np.sign(X[n][0] - 1.4) == np.sign(X[anc][0] - 1.4)}
            new = min(pool, key=lambda n: np.linalg.norm(X[n] - X[anc]))
            sb = sb.replace(f'<delem>{anc},{old}</delem>', f'<delem>{anc},{new}</delem>')
            q = arc_point(O, X[new])
            anchors[anc] = q
            targets.add(new)
            blk, k = re.subn(r'\s*<linear_constraint>\s*<node id="%d" bc="[xyz]">1</node>.*?</linear_constraint>' % anc, '',
                             blk, flags=re.S)
            assert k == 3, (anc, k)
            rep3.append(f'{nm} {anc}->{old} => {anc}->{new} (curve {half} half; anchor moved {np.linalg.norm(q - X[anc]):.2f} mm '
                        f'to the outer arc; length {np.linalg.norm(X[anc] - X[old]):.2f} -> {np.linalg.norm(q - X[new]):.2f} mm)')
        t = t[:i] + sb + t[j:]
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    t = t[:i0] + blk + t[i1:]
    for n, q in anchors.items():
        pat = f'<node id="{n}">[^<]*</node>'
        assert len(re.findall(pat, t)) == 1
        t = re.sub(pat, f'<node id="{n}">{q[0]:.8g},{q[1]:.8g},{q[2]:.8g}</node>', t)
    # 4. new springs from the distal end face's free outer nodes (outer layer: outer_end, above)
    free = sorted(n for n in outer_end if n not in wall and n not in tied and n not in targets)
    nid = max(X) + 1
    new_nodes, delems, rep4 = [], [], []
    for n in free:
        q = arc_point(O, X[n])
        new_nodes.append((nid, q))
        delems.append((n, nid))
        half = 'apex' if n in mem[CURVES[0]] and n in mem[CURVES[1]] else ('AVW' if n in mem[CURVES[0]] else 'PVW')
        rep4.append(f'{n} ({"left" if X[n][0] < 1.4 else "right"}, {half}, r {rad[n][0]:.2f} / outer {rad[n][1]:.2f}) -> '
                    f'{nid}: {np.linalg.norm(q - X[n]):.2f} mm')
        nid += 1
    j = t.rindex('</Nodes>') + len('</Nodes>')
    t = (t[:j] + '\n\t\t<Nodes name="curve_end_PM_anchor_nodes">\n'
         + ''.join(f'\t\t\t<node id="{k}">{q[0]:.8g},{q[1]:.8g},{q[2]:.8g}</node>\n' for k, q in new_nodes) + '\t\t</Nodes>' + t[j:])
    i, j = block(t, 'DiscreteSet', 'PM_conn')
    t = (t[:j] + '\n\t\t<DiscreteSet name="curve_end_PM_conn">\n'
         + ''.join(f'\t\t\t<delem>{a},{b}</delem>\n' for a, b in delems) + '\t\t</DiscreteSet>' + t[j:])
    m = re.search(r'<discrete dmat="(\d+)" discrete_set="PM_conn" />', t)
    t = t[:m.end()] + f'\n\t\t<discrete dmat="{m.group(1)}" discrete_set="curve_end_PM_conn" />' + t[m.end():]
    fixed = sorted(anchors) + [k for k, _ in new_nodes]
    j = t.index('<NodeSet name="BC-PM-outer-arc">')
    j = t.index('</NodeSet>', j) + len('</NodeSet>')
    t = t[:j] + '\n\t\t<NodeSet name="curve_end_PM_anchors">' + ','.join(map(str, fixed)) + '</NodeSet>' + t[j:]
    j = t.index('<bc name="BC-PM-outer-arc"')
    j = t.index('</bc>', j) + len('</bc>')
    t = (t[:j] + '\n\t\t<bc name="curve_end_PM_anchors" node_set="curve_end_PM_anchors" type="zero displacement">'
         '\n\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n\t\t\t<z_dof>1</z_dof>\n\t\t</bc>' + t[j:])
    msg = ('NOT IN SOURCE: the 4 PM_avw_bottom springs (the PM-LA law at full) between the curves\' distal ends and the PM\'s '
           'outer arc, their PM ends fixed there (bc curve_end_PM_anchors), their PM_anchor_ties entries dropped '
           f'({n0} -> {n0 - 12} linear constraints): ' + '; '.join(rep3)
           + f'\nNOT IN SOURCE: {len(delems)} new springs (set curve_end_PM_conn, material PM_conn_mat as PM_conn) from the '
           'outer-surface nodes of the curves\' distal end faces not on the AVW / PVW and without a spring, to new nodes at '
           'the nearest point on the PM\'s outer arc, fixed (bc curve_end_PM_anchors): ' + '; '.join(rep4))
    print(msg)
    write(ENDSPR, t, f'base: claude_diag/runs/{OUTER}/{OUTER}.feb (tools/build_batch166.py)\n{msg}\n')


def step56():
    p = os.path.join(RUNS, ENDSPR, ENDSPR + '.feb')
    t, msg5 = add_self_contact(open(p, encoding='utf-8').read())
    t, msg6 = wall_damping(t)
    write(SELF, t, f'base: claude_diag/runs/{ENDSPR}/{ENDSPR}.feb (tools/build_batch166.py)\n{msg5}\n{msg6}\n')
    print(' '.join(build(166, 'wrapfull_endspr_walls035', {'walls': 0.35}, ['healthy', 'P1', 'P2'], SELF)))


if __name__ == '__main__':
    step1()
    step2()
    step34()
    step56()
