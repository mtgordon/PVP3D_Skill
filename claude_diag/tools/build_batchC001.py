"""Batch C001 (2026-10-24, Stephen's turn; Mark's e-mail: "take the levator ani from your model and use it in all
future models"). Puts Stephen's L91_v18 levator ("LA package") onto any run of this line. Same mesh: the levator's
1170 shells / 1233 nodes and the 62 arcus-chain nodes 23677-23738 are node-for-node identical to L91_v18 (checked, coords
within 1e-6 mm), so nothing is remeshed. Three parts, applied together as ONE change (Stephen's levator):

 1. Arcus anchored to the pelvic sidewall: the 62 ATLA / posterior-arcus chain nodes fixed (L91 v2). Before, the chains
    hang on their end nodes and the ATLA / P-arcus springs only, so the levator's sides can swing out ("ballooning").
    NOT IN SOURCE (the source holds the chains at their two ends only). Side effect: the 26 Parcus_conn ends on the
    chains become fixed points.
 2. Levator healthy law = the L91 "100 % Yamada" Yeoh (c1 0.04694256, c2 0.06566996, k 4) = 2x this line's healthy law
    (build_batch107.LA100, c1 0.02347). The support case's levator % is KEPT and applied to the new law (MAXl_ant1:
    40 % -> c1 0.01878). NOT FAITHFUL (2x the source's LA_Yamada50% fit). The material keeps its name
    (LA_Yamada50pct_Yeoh) so the other tools still find it; build_batch424 asserts the OLD law, so switch cases BEFORE
    this script.
 3. Wall-levator contacts (PVW_LA, walls_LA): penalty 0.5 -> 2, search_tol 0.01 -> 0.1 (L91 v9-v12; without them the
    L242 port stalled at t 0.66). Surfaces unchanged (walls_LA already holds the AVW / cervix / curve faces; every
    canal node starts >= 3.5 mm from the LA mid-surface, so nothing starts inside the levator).

usage: py -3.10 build_batchC001.py SRC_RUN CASE NEW_RUN   (CASE = the source run's support case)
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch244 import CASES  # noqa: E402
from build_batch107 import LA100  # noqa: E402
from build_batch162 import write  # noqa: E402
from paths import RUNS_DIR  # noqa: E402

LA_L91 = {'c1': 0.04694256, 'c2': 0.06566996, 'k': 4.0}
ARCUS = list(range(23677, 23739))
BC = ('\t\t<bc name="BC-ATLA-Arcus-sidewall" node_set="BC-ATLA-Arcus-sidewall" type="zero displacement">\n'
      '\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n\t\t\t<z_dof>1</z_dof>\n\t\t</bc>\n')


def convert(t, case):
    to = case
    log = []
    assert 'BC-ATLA-Arcus-sidewall' not in t, 'already has the L91 levator'
    # 1. arcus anchoring
    for i in ARCUS:
        assert f'<node id="{i}">' in t, i
    rows = [', '.join(map(str, ARCUS[k:k + 8])) for k in range(0, len(ARCUS), 8)]
    ns = '\t\t<NodeSet name="BC-ATLA-Arcus-sidewall">\n' + ''.join(f'\t\t\t{r},\n' for r in rows)[:-2] + '\n\t\t</NodeSet>\n'
    m = re.search(r'\n(\t\t<NodeSet name="BC-ATLA-L">)', t)
    t = t[:m.start(1)] + ns + t[m.start(1):]
    m = re.search(r'\n(\t\t<bc name="BC-ATLA-L")', t)
    t = t[:m.start(1)] + BC + t[m.start(1):]
    log.append('arcus anchored: 62 ATLA / P-arcus chain nodes 23677-23738 fixed (BC-ATLA-Arcus-sidewall) [NOT IN SOURCE]')
    # 2. levator law (keeps the case's levator %)
    pa, pb = CASES[case][0], CASES[to][0]
    m = re.search(r'<material id="\d+" name="LA_Yamada50pct_Yeoh" type="Yeoh">.*?</material>', t, re.S)
    blk = m.group(0)
    for k, v in LA100.items():
        old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), blk).group(1))
        assert abs(old - v * pa / 100) < 1e-6 * max(1, v), (k, old, 'not case', case)
        blk = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{LA_L91[k] * pb / 100:.7g}</{k}>', blk)
    t = t.replace(m.group(0), blk)
    log.append(f'levator law -> L91 v18 Yeoh (c1 {LA_L91["c1"]:g}, c2 {LA_L91["c2"]:g}, k {LA_L91["k"]:g}) = 2x LA100, '
               f'at the case\'s {pb} % (was {pa} % of LA100) [NOT FAITHFUL]')
    # 3. wall-levator contacts
    for nm in ('PVW_LA', 'walls_LA'):
        m = re.search(r'<contact name="%s" .*?</contact>' % nm, t, re.S)
        blk = m.group(0)
        assert '<penalty>0.5</penalty>' in blk and '<search_tol>0.01</search_tol>' in blk, nm
        new = blk.replace('<penalty>0.5</penalty>', '<penalty>2</penalty>').replace(
            '<search_tol>0.01</search_tol>', '<search_tol>0.1</search_tol>')
        t = t.replace(blk, new)
    log.append('contacts PVW_LA + walls_LA: penalty 0.5 -> 2, search_tol 0.01 -> 0.1 (as L91 v12) [NOT IN SOURCE]')
    return t, log


if __name__ == '__main__':
    src, case, new = sys.argv[1:4]
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8', newline='').read()
    t, log = convert(t, case)
    print(new); print('\n'.join(log))
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batchC001.py, case {case})\n'
                  'ONE change: Stephen\'s L91_v18 levator (Mark 10-24: "use it in all future models"):\n- '
                  + '\n- '.join(log) + '\n')
