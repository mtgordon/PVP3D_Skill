"""Batch 293 (2026-10-16; the user, on L288c_springpatch_edge_MAX: "some elements that end up very large/oddly shaped" (the curve
tip hexes 27822-27830, 28678, 28687: V/V0 up to 7.5 vs 3.3 in L283) and the two curve_end_PM_conn springs E29213 / E29216
(anchors 29337 / 29362 on nodes 25898 / 27586) pulling out through the tip's side; "try both" as "two files").
  tipk   NOT IN SOURCE: curve hexes (canal_wrap_hex_avw / _pvw / _pvw_fine) with centre at rest z < ZCUT move to their own domains
         *_tip on Curve_soft_tip = Curve_soft with k x KF (volume stiffness only; c1 / c2 unchanged); same mass damping.
  wide   the two curve-end springs' tip-edge patches (anchor 29337 on 25898, 29362 on 27586; node 1/2 + edge neighbours 1/4)
         widened to 5 nodes along the tip edge: node 3/9, edge neighbours 2/9, next ones 1/9.
usage: py -3.10 build_batch293.py SRC_RUN NEW_RUN tipk [--zcut -38] [--kf 5]   |   SRC_RUN NEW_RUN wide"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch278 import lc  # noqa: E402
from build_batch162 import write  # noqa: E402

WIDE = {29337: [(25898, 3 / 9), (25897, 2 / 9), (25899, 2 / 9), (25896, 1 / 9), (25900, 1 / 9)],
        29362: [(27586, 3 / 9), (27585, 2 / 9), (27587, 2 / 9), (27584, 1 / 9), (27588, 1 / 9)]}


def tipk(t, X, E, zcut, kf):
    m = re.search(r'\t\t<material id="(\d+)" name="Curve_soft" type="Yeoh">.*?</material>\n', t, re.S)
    k0 = float(re.search(r'<k>([^<]+)</k>', m.group(0)).group(1))
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    mat = (m.group(0).replace(f'id="{m.group(1)}"', f'id="{mid}"').replace('name="Curve_soft"', 'name="Curve_soft_tip"')
           .replace(re.search(r'<k>[^<]+</k>', m.group(0)).group(0), f'<k>{k0 * kf:g}</k>'))
    t = t.replace(m.group(0), m.group(0) + mat, 1)
    n = {}
    for dom in ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine'):
        ids = [k for k, (nm, c) in E.items() if nm == dom and np.mean([X[v][2] for v in c]) < zcut]
        a = t.index(f'<Elements type="hex8" name="{dom}">'); b = t.index('</Elements>', a)
        blk = t[a:b]; rows = []
        for k in ids:
            row = f'\t\t\t<elem id="{k}">' + ','.join(map(str, E[k][1])) + '</elem>\n'
            assert row in blk, k
            blk = blk.replace(row, '', 1); rows.append(row)
        t = t[:a] + blk + t[b:]
        a = t.index(f'<Elements type="hex8" name="{dom}">'); b = t.index('</Elements>', a) + len('</Elements>\n')
        t = t[:b] + f'\t\t<Elements type="hex8" name="{dom}_tip">\n' + ''.join(rows) + '\t\t</Elements>\n' + t[b:]
        t = t.replace(f'\t\t<SolidDomain name="{dom}" mat="Curve_soft" />\n',
                      f'\t\t<SolidDomain name="{dom}" mat="Curve_soft" />\n\t\t<SolidDomain name="{dom}_tip" mat="Curve_soft_tip" />\n', 1)
        bl = re.search(rf'\t\t<body_load name="wall_damping_{dom}" .*?</body_load>\n', t, re.S).group(0)
        t = t.replace(bl, bl + bl.replace(dom, f'{dom}_tip'), 1)
        n[dom] = len(ids)
    return t, (f'NOT IN SOURCE: curve hexes with centre at rest z < {zcut} ({n}) moved to *_tip domains on Curve_soft_tip = '
               f'Curve_soft with k x{kf:g} ({k0:g} -> {k0 * kf:g}; c1 / c2 unchanged), same mass damping')


def wide(t):
    blk = re.search(r'\t\t<constraint name="spring_patch_ties".*?</constraint>\n', t, re.S).group(0)
    new = blk
    for a, p in WIDE.items():
        rows = [r for r in re.findall(r'\t\t\t<linear_constraint>\n.*?\t\t\t</linear_constraint>\n', blk, re.S)
                if re.search(rf'<node id="{a}" bc="[xyz]">1</node>', r)]
        assert len(rows) == 3, (a, len(rows))
        for r in rows:
            new = new.replace(r, '', 1)
        new = new.replace('\t\t</constraint>\n', lc(a, p) + '\t\t</constraint>\n')
    t = t.replace(blk, new, 1)
    return t, ('NOT IN SOURCE: the two curve_end_PM_conn tip-edge patches (anchor 29337 on 25898, 29362 on 27586) widened to 5 '
               'nodes along the tip edge (3/9, 2/9, 2/9, 1/9, 1/9) instead of 3 (1/2, 1/4, 1/4)')


def main(src, new, mode, zcut=-38.0, kf=5.0):
    t, X, E = load(src)
    t, msg = tipk(t, X, E, zcut, kf) if mode == 'tipk' else wide(t)
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch293.py {mode})\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], a[2], g('--zcut', -38.0), g('--kf', 5.0))
