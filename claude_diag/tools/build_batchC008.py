"""Batch C008 (2026-10-25, Stephen's long session; Mark's corner idea 5: "the strap edge may drag the curve tip into
the fold: spread / soften the strap's tip-side attachment"). curve_cross_diag.py on C004: the deepest curve self-crossings
(1.1-1.4 mm) sit at the curve tip right next to the 5 nodes the strap (PM_avw_bottom_*_sheet_solid) shares with the curve
(left 25896-25900, z -38..-41). Here the strap's penta6 elements whose centroid lies within RADIUS mm of those shared nodes
move to twin domains <strap>_tip with the strap material x SCALE (Mooney-Rivlin c1 and fibre E; k kept), so the strap's
pull reaches the curve through a softer, more compliant zone instead of a stiff edge. The strap's mat_axis ElementData is
split to match (lid order kept). NOT IN SOURCE (the strap itself is ours).
usage: py -3.10 build_batchC008.py SRC_RUN NEW_RUN [--radius 3] [--scale 0.25]"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from paths import RUNS_DIR  # noqa: E402

STRAPS = ('PM_avw_bottom_left_sheet_solid', 'PM_avw_bottom_right_sheet_solid')
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine',
          'canal_wrap_hex_avw_band', 'canal_wrap_hex_pvw_band', 'canal_wrap_hex_pvw_fine_band')


def block(t, name):
    return re.search(r'(\t\t<Elements type="(\w+)" name="%s">\n)(.*?)(\t\t</Elements>\n)' % name, t, re.S)


def convert(t, radius, scale):
    X = {int(i): np.array([float(v) for v in xyz.split(',')])
         for i, xyz in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    curve = set()
    for c in CURVES:
        m = block(t, c)
        if m:
            curve |= {int(w) for v in re.findall(r'<elem id="\d+">([^<]+)<', m.group(3)) for w in v.split(',')}
    log = []
    ids = [int(v) for v in re.findall(r'<material id="(\d+)"', t)]
    for s in STRAPS:
        m = block(t, s)
        lines = m.group(3).splitlines(keepends=True)
        conn = [[int(v) for v in re.search(r'<elem id="\d+">([^<]+)<', ln).group(1).split(',')] for ln in lines]
        shared = sorted({n for c in conn for n in c} & curve)
        S = np.array([X[n] for n in shared])
        tip = [k for k, c in enumerate(conn)
               if np.min(np.linalg.norm(S - np.mean([X[n] for n in c], 0), axis=1)) < radius]
        keep = [k for k in range(len(lines)) if k not in set(tip)]
        new = (m.group(1) + ''.join(lines[k] for k in keep) + m.group(4)
               + f'\t\t<Elements type="{m.group(2)}" name="{s}_tip">\n' + ''.join(lines[k] for k in tip) + m.group(4))
        t = t.replace(m.group(0), new)
        # material
        mm = re.search(r'<SolidDomain name="%s" mat="([^"]+)" />' % s, t)
        mat = mm.group(1)
        mb = re.search(r'(\t\t<material id="(\d+)" name="%s" type="[^"]+">.*?</material>\n)' % re.escape(mat), t, re.S).group(1)
        ids.append(max(ids) + 1)
        nb = re.sub(r'id="\d+"', f'id="{ids[-1]}"', mb, count=1).replace(f'name="{mat}"', f'name="{mat}_tip"')
        nb = re.sub(r'<c1>([^<]+)</c1>', lambda q: f'<c1>{float(q.group(1)) * scale:.7g}</c1>', nb)
        nb = re.sub(r'<E>([^<]+)</E>', lambda q: f'<E>{float(q.group(1)) * scale:.7g}</E>', nb)
        i = t.index('\t</Material>')
        t = t[:i] + nb + t[i:]
        t = t.replace(mm.group(0), mm.group(0) + f'\n\t\t<SolidDomain name="{s}_tip" mat="{mat}_tip" />')
        # mat_axis split (lid = position in the element set)
        em = re.search(r'(\t\t<ElementData type="mat_axis" elem_set="%s">\n)(.*?)(\t\t</ElementData>\n)' % s, t, re.S)
        rows = em.group(2).splitlines(keepends=True)
        assert len(rows) == len(lines), (s, len(rows), len(lines))
        relid = lambda r, k: re.sub(r'lid="\d+"', f'lid="{k}"', r)  # noqa: E731
        nd = (em.group(1) + ''.join(relid(rows[k], j + 1) for j, k in enumerate(keep)) + em.group(3)
              + f'\t\t<ElementData type="mat_axis" elem_set="{s}_tip">\n'
              + ''.join(relid(rows[k], j + 1) for j, k in enumerate(tip)) + em.group(3))
        t = t.replace(em.group(0), nd)
        log.append(f'{s}: {len(tip)} of {len(lines)} strap elements within {radius:g} mm of the {len(shared)} curve-shared '
                   f'nodes -> {s}_tip, material {mat}_tip = {mat} with c1 and fibre E x{scale:g}')
    for s in STRAPS:   # other references to the strap element sets (damping, outputs) would miss the tip elements
        refs = [r for r in re.findall(r'elem_set="%s"' % s, t)]
        assert len(refs) == 1, (s, 'referenced elsewhere', len(refs))
    return t, log


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('new')
    ap.add_argument('--radius', type=float, default=3); ap.add_argument('--scale', type=float, default=0.25)
    a = ap.parse_args()
    t = open(os.path.join(RUNS_DIR, a.src, a.src + '.feb'), encoding='utf-8', newline='').read()
    t, log = convert(t, a.radius, a.scale)
    print(a.new); print('\n'.join(log))
    write(a.new, t, f'base: claude_diag/runs/{a.src}/{a.src}.feb (tools/build_batchC008.py)\nONE change (Mark\'s corner idea 5: '
                    'soften the strap where it pulls on the curve tip) [NOT IN SOURCE]:\n- ' + '\n- '.join(log) + '\n')
