"""Batch 161 (2026-10-03, the user: "What if the edge nodes were constrained to not move in the x direction? Could the springs be
removed then?" ... "yes, run both variants").
On C3 on the fast copy (L157_C3_parcus042_<case>: L149_seamspr_pm_outer_fast + healthy Parcus x0.42, then the case):
  xfix         the canal's seam edge nodes (every node in constraint canal_seam_springs: the 116 AVW / cervix edge nodes and the
               PVW edge nodes they are tied to) held in x (bc seam_edge_xfix, zero displacement x only); the seam springs kept.
               NOT IN SOURCE (the source has no such boundary condition; the edges move 2-6 mm in x in C3).
  xfix_nospr   the same, and the seam springs (constraint canal_seam_springs) removed: the seam held only by the lumen contact
               SlidingElastic1 and the x constraint.
  L161_C3_xfix_<case>, L161_C3_xfix_nospr_<case>    (case = healthy, P1, P2)
usage: py -3.10 build_batch161.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')

if __name__ == '__main__':
    for case in ('healthy', 'P1', 'P2'):
        src = f'L157_C3_parcus042_{case}'
        p = os.path.join(RUNS, src, src + '.feb')
        t0 = open(p, encoding='ISO-8859-1').read()
        i = t0.index('<constraint name="canal_seam_springs"')
        j = t0.index('</constraint>', i) + len('</constraint>')
        nodes = sorted(set(int(n) for n in re.findall(r'<node id="(\d+)" bc="x">', t0[i:j])))
        # the x constraint
        k = t0.index('<NodeSet name="BC-PM-outer-arc">')
        k = t0.index('</NodeSet>', k) + len('</NodeSet>')
        t1 = t0[:k] + '\n\t\t<NodeSet name="seam_edge_nodes">' + ','.join(map(str, nodes)) + '</NodeSet>' + t0[k:]
        k = t1.index('</Boundary>')
        t1 = (t1[:k] + '\t<bc name="seam_edge_xfix" node_set="seam_edge_nodes" type="zero displacement">\n\t\t\t<x_dof>1</x_dof>'
              '\n\t\t\t<y_dof>0</y_dof>\n\t\t\t<z_dof>0</z_dof>\n\t\t</bc>\n\t' + t1[k:])
        # without the seam springs
        i = t1.index('<constraint name="canal_seam_springs"')
        j = t1.index('</constraint>', i) + len('</constraint>')
        i0 = t1.rindex('\n', 0, i) + 1
        t2 = t1[:i0] + t1[j:].lstrip(' \t').lstrip('\n')
        assert 'canal_seam_springs' not in t2
        prev = open(p + '.changes.txt', encoding='utf-8', errors='replace').read()
        msgx = (f'NOT IN SOURCE: the canal\'s {len(nodes)} seam edge nodes (every node in canal_seam_springs: AVW / cervix edge '
                f'nodes and their PVW partners) held in x (bc seam_edge_xfix, zero displacement x)')
        for name, txt, msg in ((f'L161_C3_xfix_{case}', t1, msgx + '; the seam springs kept'),
                               (f'L161_C3_xfix_nospr_{case}', t2, msgx + '; the seam springs (canal_seam_springs, 348 linear '
                                'constraints, penalty 0.01) removed')):
            d = os.path.join(RUNS, name)
            assert not os.path.exists(d), f'{name} exists; not overwriting'
            os.makedirs(d)
            open(os.path.join(d, name + '.feb'), 'w', encoding='ISO-8859-1', newline='').write(txt)
            with open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
                fo.write(f'{name} (tools/build_batch161.py) from {p}:\n{prev}\n{msg}\n')
            print('wrote', name, len(nodes), 'nodes')
