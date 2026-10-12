"""Batch 163 (2026-10-05, the user, after L162_thin275_wrapend_outer_avwbcurve_full_P1 (the 4 springs barely stretched:
their PM ends, tied to the deformable PM shell, moved ~8 mm with the walls): "fix the spring ends where they are, like the
source. Also, it looks like the attachment to the AVW with a negative x value isn't attached to the node on the end of the AVW.
Please move it to there, fix the end of the springs, and then rerun P1").
From L162_thin275_wrapend_outer_avwbcurve_full (the 4 PM_avw_bottom springs on the curves, the PM-LA law at full):
  * the left upper spring (PM_avw_bottom_left_conn 10607) moved from curve node 25868 (one row in from the curve's end) to
    25897, the outer-surface node on the end face of the left curve's AVW half (as the right one, 27553, already is);
  * the 4 PM ends (10607-10610 = the source's assembly nodes 34-37, Set-PM-AVW-bottom-origins) fixed in x, y, z as the
    source's BC-RP-PM-AVW-bottom-origins (PINNED), instead of tied to the PM shell: their PM_anchor_ties entries dropped.
  L163_thin275_wrapend_avwbend_full_fixed      the model (healthy)
  L163_thin275_wrapend_avwbend_full_fixed_P1   + the paper's P1 (build_batch145 --cases P1)
usage: py -3.10 build_batch163.py"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch145 import build  # noqa: E402
from febmodel import Feb  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L162_thin275_wrapend_outer_avwbcurve_full'
NAME = 'L163_thin275_wrapend_avwbend_full_fixed'
ANCHORS = (10607, 10608, 10609, 10610)
MOVE = (10607, 25868, 25897)


def main():
    p = os.path.join(RUNS, SRC, SRC + '.feb')
    t = open(p, encoding='utf-8').read()
    X = {n: np.array(v, float) for n, v in Feb(p).nodes.items()}
    a, old, new = MOVE
    assert t.count(f'<delem>{a},{old}</delem>') == 1
    t = t.replace(f'<delem>{a},{old}</delem>', f'<delem>{a},{new}</delem>')
    msg1 = (f'NOT IN SOURCE: PM_avw_bottom_left_conn {a}->{old} ({np.linalg.norm(X[a] - X[old]):.2f} mm, one row in from the '
            f'curve\'s end) => {a}->{new} ({np.linalg.norm(X[a] - X[new]):.2f} mm; node at '
            f'{", ".join(f"{v:.2f}" for v in X[new])}: the outer-surface node on the end face of the left curve\'s AVW half, as '
            f'the right spring\'s 27553)')
    i0 = t.index('<constraint name="PM_anchor_ties"')
    i1 = t.index('</constraint>', i0)
    blk = t[i0:i1]
    n0 = blk.count('<linear_constraint>')
    for n in ANCHORS:
        blk, k = re.subn(r'\s*<linear_constraint>\s*<node id="%d" bc="[xyz]">1</node>.*?</linear_constraint>' % n, '', blk,
                         flags=re.S)
        assert k == 3, (n, k)
    t = t[:i0] + blk + t[i1:]
    assert 'BC-RP-PM-AVW-bottom-origins' not in t
    j = t.index('<NodeSet name="BC-PM-outer-arc">')
    j = t.index('</NodeSet>', j) + len('</NodeSet>')
    t = t[:j] + '\n\t\t<NodeSet name="BC-RP-PM-AVW-bottom-origins">' + ','.join(map(str, ANCHORS)) + '</NodeSet>' + t[j:]
    j = t.index('<bc name="BC-PM-outer-arc"')
    j = t.index('</bc>', j) + len('</bc>')
    t = (t[:j] + '\n\t\t<bc name="BC-RP-PM-AVW-bottom-origins" node_set="BC-RP-PM-AVW-bottom-origins" type="zero displacement">'
         '\n\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n\t\t\t<z_dof>1</z_dof>\n\t\t</bc>' + t[j:])
    msg2 = (f'as the source: the 4 PM_avw_bottom PM ends {", ".join(map(str, ANCHORS))} (source assembly nodes 34-37, '
            f'Set-PM-AVW-bottom-origins) fixed in x, y, z (bc BC-RP-PM-AVW-bottom-origins, zero displacement; the source\'s '
            f'PINNED) where they are, instead of tied to the PM shell: their PM_anchor_ties entries dropped ({n0} -> '
            f'{n0 - 3 * len(ANCHORS)} linear constraints)')
    print(msg1)
    print(msg2)
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch163.py)\n{msg1}\n{msg2}\n')
    print(' '.join(build(163, NAME[len('L163_'):], {}, ['P1'], NAME)))


if __name__ == '__main__':
    main()
