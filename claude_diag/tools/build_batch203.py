"""Batch 203 (2026-10-07, the user: "Test out your idea for the tube line"). On the tube the AVW and PVW descend together
(posterior x0.4: Ba +6.9 with Bp +9.5; a softer PVW-side curve dropped both). On the band line the gap opened only once
the seam was contact-only. The tube's equivalent: soften the join between the curves and the PVW. One change on the tube
default L202_curvepress_apexclosed_P1 (closed apex, whole-curve pressure, curve damping C 1000, equal PM-arch springs), P1:
  L203_apexclosed_pvwsoftjoin_P1   the 340 curve nodes tied to the PVW's columns (canal_wrap_ends rows whose wall nodes are
                                   PVW: 1020 of 1692) move to their own constraint canal_wrap_pvw_soft, penalty 10 -> 0.05
                                   N/mm per dof (about the wall's own nodal stiffness), maxaug 10 -> 0 (an elastic
                                   zero-length spring, not an exact tie). NOT IN SOURCE. The 200 nodes the PVW half SHARES
                                   with the PVW stay rigid (not duplicated here).
usage: py -3.10 build_batch203.py [PENALTY]"""
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L202_curvepress_apexclosed_P1'

if __name__ == '__main__':
    pen = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
    name = 'L203_apexclosed_pvwsoftjoin_P1' if pen == 0.05 else f'L203_apexclosed_pvwsoftjoin{pen:g}_P1'
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    pvw = {int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet64' for el in e for v in el.text.split(',')}
    a = t.index('<constraint name="canal_wrap_ends"')
    b = t.index('</constraint>', a)
    blk = t[a:b]
    keep, soft = [], []
    for m in re.finditer(r'\t*<linear_constraint>.*?</linear_constraint>\n?', blk, re.S):
        ids = [int(v) for v in re.findall(r'node id="(\d+)"', m.group(0))]
        (soft if all(i in pvw for i in ids[1:]) else keep).append(m.group(0))
    head = blk[:blk.index('<linear_constraint>')].rstrip('\t')
    t = t[:a] + head + ''.join(keep) + '\t\t' + t[b:]
    i = t.index('</constraint>', t.index('<constraint name="canal_wrap_ends"')) + len('</constraint>\n')
    t = t[:i] + (f'\t\t<constraint name="canal_wrap_pvw_soft" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n'
                 f'\t\t\t<penalty>{pen:g}</penalty>\n\t\t\t<maxaug>0</maxaug>\n' + ''.join(soft) + '\t\t</constraint>\n') + t[i:]
    write(name, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch203.py)\n'
                   f'NOT IN SOURCE: the curves\' ties to the PVW columns ({len(soft)} of {len(soft) + len(keep)} canal_wrap_ends rows) '
                   f'moved to canal_wrap_pvw_soft, penalty {pen:g} N/mm per dof, maxaug 0 (elastic, not exact); the 200 nodes '
                   f'shared by the PVW half and the PVW stay rigid\n')
    print(len(soft), 'soft rows,', len(keep), 'kept')
