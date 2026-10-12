"""Batch 204 (2026-10-07; the user's plan: "duplicate the 200 shared curve-PVW nodes for a fully elastic join").
L203 (the curves' 340 TIED nodes on the PVW made elastic, penalty 0.05, maxaug 0) matched the tube default L202 within
0.2 mm to t 0.73: the 200 nodes the curves SHARE with the PVW (_PickedSet64) still join them rigidly. One change on
L203_apexclosed_pvwsoftjoin_P1, P1:
  L204_apexclosed_pvwsoftjoin_full_P1   every node shared by a curve element (canal_wrap_hex_avw / _pvw / *_apex) and a PVW
                                        element is split: the curves get a copy (same position), every surface face that is
                                        a curve-element face takes the copy, and the copy follows the PVW node through the
                                        same elastic tie (canal_wrap_pvw_soft: penalty 0.05 N/mm per dof, maxaug 0).
                                        Springs, BCs, constraints and the PVW's own faces keep the original node. NOT IN SOURCE.
usage: py -3.10 build_batch204.py"""
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L203_apexclosed_pvwsoftjoin_P1'
NAME = 'L204_apexclosed_pvwsoftjoin_full_P1'
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]

if __name__ == '__main__':
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    E = {e.get('name'): [[int(v) for v in el.text.split(',')] for el in e] for e in mesh.findall('Elements')}
    X = {int(n.get('id')): n.text for b in mesh.findall('Nodes') for n in b}
    curve_doms = [d for d in E if d.startswith('canal_wrap_hex')]
    cn = {n for d in curve_doms for el in E[d] for n in el}
    pn = {n for el in E['_PickedSet64'] for n in el}
    S = sorted(cn & pn)
    nid = max(X) + 1
    copy = {s: nid + i for i, s in enumerate(S)}
    cfaces = {frozenset(el[i] for i in fc) for d in curve_doms for el in E[d] for fc in HF}
    # 1. curve elements take the copies
    for d in curve_doms:
        a = t.index(f'<Elements type="hex8" name="{d}">')
        b = t.index('</Elements>', a)
        blk = re.sub(r'(<elem id="\d+">)([^<]+)(</elem>)',
                     lambda m: m.group(1) + ','.join(str(copy.get(int(v), int(v))) for v in m.group(2).split(',')) + m.group(3),
                     t[a:b])
        t = t[:a] + blk + t[b:]
    # 2. surface faces that are curve-element faces take the copies
    nsw = 0

    def face(m):
        global nsw
        ids = [int(v) for v in m.group(2).split(',')]
        if frozenset(ids) in cfaces and any(i in copy for i in ids):
            nsw += 1
            ids = [copy.get(i, i) for i in ids]
        return m.group(1) + ','.join(map(str, ids)) + m.group(3)
    a = t.index('<Surface ')
    b = t.rindex('</Surface>') + len('</Surface>')
    t = t[:a] + re.sub(r'(<(?:quad4|tri3) id="\d+">)([^<]+)(</(?:quad4|tri3)>)', face, t[a:b]) + t[b:]
    # 3. the copies' nodes, after the last <Nodes> block (ids stay increasing)
    blk = '\t\t<Nodes name="curve_pvw_split_nodes">\n' + ''.join(
        f'\t\t\t<node id="{copy[s]}">{X[s]}</node>\n' for s in S) + '\t\t</Nodes>\n'
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + blk + t[i:]
    # 4. elastic ties copy <-> original, in canal_wrap_pvw_soft
    ties = ''.join(f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{copy[s]}" bc="{c}">1</node>\n\t\t\t\t<node id="{s}" bc="{c}">-1</node>\n'
                   f'\t\t\t</linear_constraint>\n' for s in S for c in 'xyz')
    a = t.index('<constraint name="canal_wrap_pvw_soft"')
    b = t.index('</constraint>', a)
    t = t[:b] + ties + '\t\t' + t[b:]
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch204.py)\n'
                   f'NOT IN SOURCE: the {len(S)} nodes shared by the curves and the PVW split (copies {nid}-{nid + len(S) - 1} '
                   f'for the curve elements and {nsw} curve surface faces), each copy tied elastically to its PVW node in '
                   f'canal_wrap_pvw_soft (penalty 0.05, maxaug 0): the whole curve-PVW join is now elastic\n')
    print(len(S), 'nodes split,', nsw, 'surface faces switched')
