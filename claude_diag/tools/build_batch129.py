"""Batch 129 (2026-09-30, the user: "also attach the PM_conn_mat to the inside arc of the pm_plane. Just shorten them as
needed to attach them to the material"). On the base L128_tube_r30s_pm (the tube + the PM as a structure):
the 26 PM_conn springs (material PM_conn_mat) run 14.0 mm from an anchor on the PM band (0.2-2.3 mm inside its outer arc,
tied to the PM facet under it by PM_anchor_ties) to an AVW node (_PickedSet347); every one crosses the PM's inner arc
(the U's inner boundary edge, build_batch113.pm_arcs) within 0.1-0.2 mm. Each spring's anchor node is moved along the
spring to that crossing (the closest point of the spring line to the inner-arc polyline), so the spring now starts on the
inner arc and is shorter (7.4-10.3 mm). The anchor's PM_anchor_ties entry is replaced by a tie to the inner-arc edge it
lies on: u_anchor - (1-u) u_k - u u_k+1 = 0 per dof (constraint PM_conn_inner_arc, penalty 100, maxaug 0, as the other
anchors). The spring law is unchanged (PM_conn_mat: force against elongation in mm, from the new rest length).
NOT IN SOURCE (the source's PM_conn connectors end at fixed points 14 mm from the AVW).
  L129_tube_r30s_pm_conninner   L128_tube_r30s_pm + the PM_conn springs attached to the PM's inner arc
usage: py -3.10 build_batch129.py"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS  # noqa: E402
from build_batch113 import pm_arcs  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L128_tube_r30s_pm'
NAME = 'L129_tube_r30s_pm_conninner'
SETS = ('PM_conn', 'PM-LA-x%stiff_mat')   # L129 itself was built with PM_conn only


def seg_seg(p0, p1, q0, q1):
    """closest points of segments p0-p1 and q0-q1: (distance, s on p, u on q)"""
    d1, d2, r = p1 - p0, q1 - q0, p0 - q0
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    c, b = d1 @ r, d1 @ d2
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-12 else 0.0
    u = (b * s + f) / e
    if u < 0:
        u, s = 0.0, np.clip(-c / a, 0, 1)
    elif u > 1:
        u, s = 1.0, np.clip((b - c) / a, 0, 1)
    return np.linalg.norm(p0 + s * d1 - (q0 + u * d2)), s, u


if __name__ == '__main__':
    src = os.path.join(RUNS, BASE, BASE + '.feb')
    d = os.path.join(RUNS, NAME)
    assert not os.path.exists(d), f'{NAME} exists; not overwriting'
    m = Model7(src)
    X = m.nodes()
    _, _, inner, _ = pm_arcs(src)
    t = open(src, encoding='utf-8').read()
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    tied = set(int(v) for v in re.findall(r'<node id="(\d+)" bc="x">1<', t[i0:i1]))
    # the 26 PM_conn springs and the 27th connector of the same source family (PM-LA-x%stiff: the middle one, which the
    # converter wrote as its own set PM-LA-x%stiff_mat); added 2026-10-01 at the user's request
    dss = [s for s in m.mesh.findall('DiscreteSet') if s.get('name') in SETS]
    moves, ties, rep = {}, [], []
    for e in [e for ds in dss for e in ds.findall('delem')]:
        a, b = (int(v) for v in e.text.split(','))
        anc, tis = (a, b) if a in tied else (b, a)
        assert anc in tied and tis not in tied
        A, T = np.array(X[anc]), np.array(X[tis])
        best = min((seg_seg(A, T, np.array(X[inner[k]]), np.array(X[inner[k + 1]])) + (k,) for k in range(len(inner) - 1)),
                   key=lambda r: r[0])
        dist, s, u, k = best
        assert dist < 0.5, (anc, dist)
        q = np.array(X[inner[k]]) * (1 - u) + np.array(X[inner[k + 1]]) * u
        moves[anc] = q
        ties.append((anc, [(inner[k], 1 - u), (inner[k + 1], u)]))
        rep.append((np.linalg.norm(T - A), np.linalg.norm(T - q), dist))
    # 1. move the anchor nodes onto the inner arc
    for n, q in moves.items():
        pat = f'<node id="{n}">[^<]*</node>'
        assert len(re.findall(pat, t)) == 1, n
        t = re.sub(pat, f'<node id="{n}">{q[0]:.8g},{q[1]:.8g},{q[2]:.8g}</node>', t)
    # 2. drop their PM_anchor_ties entries
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    blk = t[i0:i1]
    n_before = blk.count('<linear_constraint>')
    for n in moves:
        blk, k = re.subn(r'\s*<linear_constraint>\s*<node id="%d" bc="[xyz]">1</node>.*?</linear_constraint>' % n, '', blk,
                         flags=re.S)
        assert k == 3, (n, k)
    t = t[:i0] + blk + t[i1:]
    # 3. tie them to the inner-arc edge
    lines = ['\t\t<constraint name="PM_conn_inner_arc" type="linear constraint">', '\t\t\t<tol>0.01</tol>',
             '\t\t\t<penalty>100</penalty>', '\t\t\t<maxaug>0</maxaug>']
    for n, ws in ties:
        for dof in 'xyz':
            lines.append('\t\t\t<linear_constraint>')
            lines.append(f'\t\t\t\t<node id="{n}" bc="{dof}">1</node>')
            for q_, w in ws:
                if w > 1e-9:
                    lines.append(f'\t\t\t\t<node id="{q_}" bc="{dof}">{-w:.6f}</node>')
            lines.append('\t\t\t</linear_constraint>')
    lines.append('\t\t</constraint>\n')
    j = t.index('</Constraints>')
    t = t[:j] + '\n'.join(lines) + '\t' + t[j:]
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t)
    L0, L1, dd = (np.array(v) for v in zip(*rep))
    msg = (f'{NAME}: {BASE} + NOT IN SOURCE: the {len(moves)} PM_conn-family springs attached to the PM\'s inner arc: each anchor node moved '
           f'along its spring to where the spring crosses the inner arc (spring-to-arc distance max {dd.max():.2f} mm), so '
           f'the springs shorten from {L0.min():.1f}-{L0.max():.1f} to {L1.min():.1f}-{L1.max():.1f} mm (PM_conn_mat '
           f'unchanged: force vs elongation from the new rest length); their PM_anchor_ties entries ({n_before} -> '
           f'{n_before - 3 * len(moves)} linear constraints) replaced by ties to the inner-arc edge (constraint '
           f'PM_conn_inner_arc, penalty 100, maxaug 0)\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
