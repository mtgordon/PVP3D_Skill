"""Batch 276 (2026-10-10, the user: "go ahead and run both of the suggested runs"). The tube main's MAX (L258_tubecand_MAX):
the curves' distal PVW-half elements next to the canal_wrap_embed ties (20 curve end nodes embedded in the PVW / PeB hexes they
lie in; build_batch143) stretch to 4-6x their volume from t 0.6. One change each:
  noembed  the canal_wrap_embed constraint removed (a diagnostic: are the ties what tears those elements? the curve ends are
           then held only by the curve and the curve-end springs; their faces are in no contact)
  spread   the tie spread: every node of the curves' PVW half (columns 3-6 of the U, all 5 rows) in the rings within 8 mm of
           the hymen (depth < 8 mm) that lies in or next to a PVW / PeB hex (|xi| <= 1.5) is embedded in it too (same
           constraint, penalty 10), so the curve end follows the PVW / PeB over a longer stretch
NOT IN SOURCE. usage: py -3.10 build_batch276.py noembed|spread"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch143 import embed  # noqa: E402
from curve_rings import rings  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC, REF = 'L258_tubecand_MAX', 'L241_curvepress_apexcontour_pmc_cx3_P1'

if __name__ == '__main__':
    mode = sys.argv[1]
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    blk = re.search(r'\t\t<constraint name="canal_wrap_embed" type="linear constraint">.*?</constraint>\n', t, re.S).group(0)
    tied = set(int(v) for v in re.findall(r'<linear_constraint>\s*<node id="(\d+)" bc="x">1</node>', blk))
    if mode == 'noembed':
        t = t.replace(blk, '')
        msg = f'DIAGNOSTIC, NOT IN SOURCE: canal_wrap_embed removed ({len(tied)} curve end nodes no longer embedded in the PVW / PeB)'
        name = 'L276_tubemain_noembed_MAX'
    else:
        mesh = ET.fromstring(t).find('Mesh')
        X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
        hosts = [([int(v) for v in el.text.split(',')], None) for e in mesh.findall('Elements')
                 if e.get('name') in ('_PickedSet64', '_PickedSet66') for el in e]
        hosts = [(c, np.array([X[n] for n in c])) for c, _ in hosts]
        HC = np.array([P.mean(0) for _, P in hosts])
        RG = rings(open(os.path.join(RUNS, REF, REF + '.feb'), encoding='utf-8').read())
        cand = set()
        for sgn, (grids, c0, B, XR) in RG.items():
            for g in grids:
                if ((XR[int(g[2, 3])] - c0) @ B)[2] < 8.0:
                    cand |= {int(v) for v in g[:, 3:].ravel()}
        new, skipped = [], 0
        for n in sorted(cand - tied):
            conn, w, mx = embed(X[n], hosts, HC)
            if mx > 1.5:
                skipped += 1
                continue
            for dof in 'xyz':
                new.append(f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{n}" bc="{dof}">1</node>\n'
                           + ''.join(f'\t\t\t\t<node id="{q}" bc="{dof}">{-wq:.6f}</node>\n' for q, wq in zip(conn, w) if abs(wq) > 1e-9)
                           + '\t\t\t</linear_constraint>\n')
        nadd = len(new) // 3
        t = t.replace(blk, blk.replace('\t\t</constraint>\n', ''.join(new) + '\t\t</constraint>\n'))
        msg = (f'NOT IN SOURCE: canal_wrap_embed spread: {nadd} more curve nodes (PVW half, rings within 8 mm of the hymen, in or '
               f'next to a PVW / PeB hex, |xi| <= 1.5) embedded in it ({len(tied)} before; {skipped} candidates too far from a host)')
        name = 'L276_tubemain_spreadembed_MAX'
    print(msg)
    write(name, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch276.py)\n{msg}\n')
