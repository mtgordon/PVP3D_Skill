"""Batch 149 (2026-10-02, the user: "I actually want the top springs that attach the AVW to the PM arch to be moved to the top
of the arch which is a boundary condition. Make a note of that for future builds (both here and going forward)").
The 27 PM_conn-family springs (the 26 PM_conn + the middle one, set PM-LA-x%stiff_mat; source law PM-LA-x%stiff) run 14.04 mm
from an AVW node to an anchor on the PM band, tied to the PM facet under it (PM_anchor_ties). Each spring line crosses the
PM's outer arc (the clamped edge, BC-PM-outer-arc) within 0.12 mm, 11.6-13.9 mm from its AVW node. Each anchor node is moved
along its spring to that crossing and fixed there (zero displacement, node set PM_conn_outer_anchors, as the clamped arc it
lies on); its PM_anchor_ties entries are dropped. The springs' law (elongation axis) is scaled by their mean new length /
14.04 mm (one factor for all 27, as build_batch129b). NOT IN SOURCE (the source's connectors end at fixed points 14 mm from
the AVW, 0.2-2.3 mm inside the outer arc).
Replaces build_batch129 + 129b (inner arc) for every new model from 2026-10-02.
  outer_arc(text, path) -> (new text, log line)     used by the builders below
  L149_seamspr_pm_outer        L144's build with PM_conn on the outer arc: L91 + seg_up 2 + PM + walls_LA / walls_PM
  L149_seamspr_pm_outer_fast   the same without walls_LA / walls_PM (as build_batch147)
usage: py -3.10 build_batch149.py"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch113 import pm_arcs  # noqa: E402
from build_batch129 import seg_seg  # noqa: E402
from febmodel import Feb  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
SETS = ('PM_conn', 'PM-LA-x%stiff_mat')
MATS = ('PM_conn_mat', 'PM-LA-x%stiff_mat')
PRE = os.path.join(SCRATCH, 'L144_seamspr_pm_allcontact_tmppre', 'L144_seamspr_pm_allcontact_tmppre.feb')


def outer_arc(t, path):
    """t: the text of the model at path (PM structure with the anchors tied by PM_anchor_ties, PM_conn on its source anchors)"""
    f = Feb(path)
    X = {n: np.array(p, float) for n, p in f.nodes.items()}
    _, outer, _, _ = pm_arcs(path)
    O = [X[n] for n in outer]
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    tied = set(int(v) for v in re.findall(r'<node id="(\d+)" bc="x">1<', t[i0:i1]))
    moves, rep = {}, []
    for nm in SETS:
        i = t.index(f'<DiscreteSet name="{nm}">')
        for a, b in re.findall(r'<delem>(\d+),(\d+)</delem>', t[i:t.index('</DiscreteSet>', i)]):
            a, b = int(a), int(b)
            anc, tis = (a, b) if a in tied else (b, a)
            assert anc in tied and tis not in tied
            A, T = X[anc], X[tis]
            E = A + (A - T) / np.linalg.norm(A - T) * 10          # the line carried 10 mm past the anchor
            dist, s, u, k = min((seg_seg(T, E, O[k], O[k + 1]) + (k,) for k in range(len(O) - 1)), key=lambda r: r[0])
            assert dist < 0.5, (anc, dist)
            q = O[k] * (1 - u) + O[k + 1] * u
            moves[anc] = q
            rep.append((np.linalg.norm(A - T), np.linalg.norm(q - T), dist))
    L0, L1, dd = (np.array(v) for v in zip(*rep))
    fac = L1.mean() / L0.mean()
    for n, q in moves.items():
        pat = f'<node id="{n}">[^<]*</node>'
        assert len(re.findall(pat, t)) == 1, n
        t = re.sub(pat, f'<node id="{n}">{q[0]:.8g},{q[1]:.8g},{q[2]:.8g}</node>', t)
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    blk = t[i0:i1]
    n_before = blk.count('<linear_constraint>')
    for n in moves:
        blk, k = re.subn(r'\s*<linear_constraint>\s*<node id="%d" bc="[xyz]">1</node>.*?</linear_constraint>' % n, '', blk,
                         flags=re.S)
        assert k == 3, (n, k)
    t = t[:i0] + blk + t[i1:]
    # the anchors fixed, as the clamped arc they now lie on
    j = t.index('<NodeSet name="BC-PM-outer-arc">')
    j = t.index('</NodeSet>', j) + len('</NodeSet>')
    t = t[:j] + '\n\t\t<NodeSet name="PM_conn_outer_anchors">' + ','.join(map(str, sorted(moves))) + '</NodeSet>' + t[j:]
    j = t.index('<bc name="BC-PM-outer-arc"')
    j = t.index('</bc>', j) + len('</bc>')
    t = (t[:j] + '\n\t\t<bc name="PM_conn_outer_anchors" node_set="PM_conn_outer_anchors" type="zero displacement">'
         '\n\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n\t\t\t<z_dof>1</z_dof>\n\t\t</bc>' + t[j:])
    for mat in MATS:
        i = t.index(f'name="{mat}" type="nonlinear spring">')
        j = t.index('</discrete_material>', i)
        blk = t[i:j]
        assert len(re.findall(r'<pt>([^,<]+),([^<]+)</pt>', blk)) == 11
        new = re.sub(r'<pt>([^,<]+),([^<]+)</pt>', lambda mo: f'<pt>{float(mo.group(1)) * fac:.6g},{mo.group(2)}</pt>', blk)
        t = t[:i] + new + t[j:]
    msg = (f'NOT IN SOURCE: the {len(moves)} PM_conn-family springs (PM_conn + PM-LA-x%stiff_mat) anchored on the PM\'s outer '
           f'arc (the clamped edge): each anchor moved along its spring to where the spring crosses the outer arc (line-to-'
           f'arc distance max {dd.max():.2f} mm) and fixed there (bc PM_conn_outer_anchors, zero displacement); springs '
           f'{L0.min():.2f}-{L0.max():.2f} -> {L1.min():.2f}-{L1.max():.2f} mm (mean {L1.mean():.2f}); their PM_anchor_ties '
           f'entries dropped ({n_before} -> {n_before - 3 * len(moves)} linear constraints); {", ".join(MATS)} elongation '
           f'axis x {fac:.4f} (mean new / old length: the same force at the same mean fractional stretch)')
    return t, msg


def strip_wall_contacts(t):
    for nm in ('walls_LA', 'walls_PM'):
        t, k = re.subn(r'[ \t]*<contact [^>]*name="%s"[^>]*>.*?</contact>\s*\n' % nm, '', t, flags=re.S)
        assert k == 1, (nm, k)
        t, k = re.subn(r'[ \t]*<SurfacePair name="%s">.*?</SurfacePair>\s*\n' % nm, '', t, flags=re.S)
        assert k == 1, (nm, k)
    for s in ('walls_LA_primary', 'walls_PM_primary', 'PM_contact'):
        t, k = re.subn(r'[ \t]*<Surface name="%s">.*?</Surface>\s*\n' % s, '', t, flags=re.S)
        assert k == 1, (s, k)
    return t


if __name__ == '__main__':
    t0 = open(PRE, encoding='ISO-8859-1').read()
    pre_log = open(PRE + '.changes.txt', encoding='utf-8').read()
    t, msg = outer_arc(t0, PRE)
    print(msg)
    for name, txt, extra in (('L149_seamspr_pm_outer', t, ''),
                             ('L149_seamspr_pm_outer_fast', strip_wall_contacts(t),
                              'the fast copy: walls_LA and walls_PM removed (contacts, SurfacePairs, surfaces walls_LA_primary / '
                              'walls_PM_primary / PM_contact), as build_batch147\n')):
        d = os.path.join(RUNS, name)
        assert not os.path.exists(d), f'{name} exists; not overwriting'
        os.makedirs(d)
        open(os.path.join(d, name + '.feb'), 'w', encoding='ISO-8859-1', newline='').write(txt)
        with open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
            fo.write(f'{name} (tools/build_batch149.py) from {PRE} (= build_batch144 before its PM_conn step):\n{pre_log}{msg}\n'
                     + extra)
        print('wrote', name)
