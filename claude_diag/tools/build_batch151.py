"""Batch 151 (2026-10-02, the user: "I don't believe the other PM Mid connections had much force at all so you might be able to
fix those or get rid of them without affecting the intermediate trials too much. Then they could be added back in at the end.")
The faster learning copy: L149_seamspr_pm_outer_fast (the outer-arc PM_conn springs, no walls_LA / walls_PM) with the 20
anchors of the near-zero PM families (PM_PeB_Left/Right_conn, PM_avw_bottom_left/right_conn: the source law PM-LA x 1e-5)
fixed where they are, as in the source, instead of tied to the PM shell: constraint PM_anchor_ties (60 linear constraints)
removed, bc PM_nearzero_anchors_fixed (zero displacement) added. The PM shell is then loaded by nothing (clamped outer arc).
For the intermediate (sensitivity) trials only; the final set goes back on L149_seamspr_pm_outer (ties restored).
  L151_seamspr_pm_outer_fast2   built, test-read, then run with its P1 / P2
usage: py -3.10 build_batch151.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE, NAME = 'L149_seamspr_pm_outer_fast', 'L151_seamspr_pm_outer_fast2'

if __name__ == '__main__':
    d = os.path.join(RUNS, NAME)
    assert not os.path.exists(d), f'{NAME} exists; not overwriting'
    t = open(os.path.join(RUNS, BASE, BASE + '.feb'), encoding='ISO-8859-1').read()
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0) + len('</constraint>')
    blk = t[i0:i1]
    anchors = sorted(set(int(v) for v in re.findall(r'<node id="(\d+)" bc="x">1<', blk)))
    nlc = blk.count('<linear_constraint>')
    assert nlc == 3 * len(anchors), (nlc, len(anchors))
    i0 = t.rindex('\n', 0, i0) + 1
    t = t[:i0] + t[i1:].lstrip(' \t').lstrip('\n')
    j = t.index('<NodeSet name="PM_conn_outer_anchors">')
    j = t.index('</NodeSet>', j) + len('</NodeSet>')
    t = t[:j] + '\n\t\t<NodeSet name="PM_nearzero_anchors">' + ','.join(map(str, anchors)) + '</NodeSet>' + t[j:]
    j = t.index('<bc name="PM_conn_outer_anchors"')
    j = t.index('</bc>', j) + len('</bc>')
    t = (t[:j] + '\n\t\t<bc name="PM_nearzero_anchors_fixed" node_set="PM_nearzero_anchors" type="zero displacement">'
         '\n\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n\t\t\t<z_dof>1</z_dof>\n\t\t</bc>' + t[j:])
    os.makedirs(d)
    open(os.path.join(d, NAME + '.feb'), 'w', encoding='ISO-8859-1', newline='').write(t)
    msg = (f'base: {os.path.join(RUNS, BASE, BASE + ".feb")} (tools/build_batch151.py)\n{NAME}: {BASE} with the {len(anchors)} '
           f'anchors of the near-zero PM families (PM-PeB, PM-AVW-bottom; PM-LA x 1e-5) fixed in place as in the source '
           f'(bc PM_nearzero_anchors_fixed) instead of tied to the PM shell (constraint PM_anchor_ties, {nlc} linear '
           f'constraints, removed); the PM shell is loaded by nothing; the learning copy for the sensitivities\n')
    open(os.path.join(d, NAME + '.feb.changes.txt'), 'w', encoding='utf-8').write(msg)
    print(msg)
