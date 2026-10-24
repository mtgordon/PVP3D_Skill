"""Batch 454 (2026-10-23; the user: "make one where the strap is replaced with connectors like the pm arch to the avw").
NOT IN SOURCE (spring layout; the strap itself was ours): the fan straps (PM_avw_bottom_left / right_sheet_solid + _wide, penta6;
build_batch337 fitted their fibre E to the 5 spring chains a side they replaced: curve_end_PM_conn x3 + PM_avw_bottom_conn x2) are
removed and replaced by single springs, like PM_conn (AVW -> top of the PM arch): one spring from each node the strap shared with
the curve / AVW (its attached edge: the curve tip -> the AVW seam node) to the strap's fixed PM-arch anchors, paired in order along
both edges (anchors picked evenly by index, no crossing). Material: a copy of PM_conn_mat's table (curve_end_PM_conn used PM_conn_mat
before the strap), <scale> = TOTAL / N per side, so the side's total equals TOTAL PM_conn springs (default 5 / 3: the 5 chains the
strap was fitted to, x1/3 as the strap since L431). Removed with the strap: its 4 element blocks and domains, their mat_axis, the
strap contacts fanstrap_curves and curve_strap (pairs + surfaces), strap_wide_ties. The strap materials stay (unused); the strap's
nodes stay (no elements; the anchors keep their zero-displacement BC).
usage: py -3.10 build_batch454.py SRC_RUN NEW_RUN [--total 1.6667]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

SIDES = ('left', 'right')
CONTACTS = ('fanstrap_curves', 'curve_strap')


def main(src, new, total):
    t, X, E = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    dom = {}
    for n, c in E.values():
        for v in c:
            dom.setdefault(v, set()).add(n)
    an = {int(v) for v in re.findall(r'\d+', re.search(r'<NodeSet name="curve_end_PM_anchors">([^<]*)<', t).group(1))}
    sets, report = '', []
    for side in SIDES:
        sn = {v for n, c in E.values() if n.startswith(f'PM_avw_bottom_{side}_sheet') for v in c}
        edge = [v for v in sn if any(d.startswith('canal_wrap') or d == '_PickedSet347' for d in dom[v])]
        anc = [v for v in sn if v in an]
        avw = [v for v in edge if '_PickedSet347' in dom[v]]
        assert len(avw) == 1, (side, avw)
        a0 = X[avw[0]]
        edge.sort(key=lambda v: np.linalg.norm(X[v] - a0))                       # from the AVW seam node to the curve tip
        start = min(anc, key=lambda v: np.linalg.norm(X[v] - a0))
        anc.sort(key=lambda v: np.linalg.norm(X[v] - X[start]))                  # from the arch end nearest the AVW
        pick = [anc[i] for i in np.round(np.linspace(0, len(anc) - 1, len(edge))).astype(int)]
        name = f'PM_avw_strap_{side}_conn'
        sets += f'\t\t<DiscreteSet name="{name}">\n' + ''.join(f'\t\t\t<delem>{e},{a}</delem>\n' for e, a in zip(edge, pick)) + '\t\t</DiscreteSet>\n'
        L = [float(np.linalg.norm(X[e] - X[a])) for e, a in zip(edge, pick)]
        report.append(f'{side}: {len(edge)} springs, edge {edge[0]}..{edge[-1]} -> anchors {len(anc)}, rest length {min(L):.1f}-{max(L):.1f} mm')
        n_side = len(edge)
    # spring material: PM_conn_mat's table, scale TOTAL / N
    dm = re.findall(r'<discrete_material id="(\d+)" name="([^"]+)"[^>]*>.*?</discrete_material>', t, re.S)
    pmc = re.search(r'<discrete_material id="\d+" name="PM_conn_mat"[^>]*>.*?</discrete_material>', t, re.S).group(0)
    nid = len(dm) + 1
    scale = total / n_side
    mat = re.sub(r'id="\d+" name="PM_conn_mat"', f'id="{nid}" name="PM_avw_strap_conn_mat"', pmc)
    mat = re.sub(r'<scale>[^<]+</scale>', f'<scale>{scale:.6g}</scale>', mat)
    last = list(re.finditer(r'\t\t<discrete_material .*?</discrete_material>\n', t, re.S))[-1]
    t = t[:last.end()] + '\t\t' + mat + '\n' + t[last.end():]
    i = t.index('\t</Discrete>')
    t = t[:i] + ''.join(f'\t\t<discrete dmat="{nid}" discrete_set="PM_avw_strap_{s}_conn" />\n' for s in SIDES) + t[i:]
    i = t.index('\t\t<DiscreteSet name=')
    t = t[:i] + sets + t[i:]
    # remove the strap
    for side in SIDES:
        for part in ('solid', 'wide'):
            d = f'PM_avw_bottom_{side}_sheet_{part}'
            t, n1 = re.subn(r'\t\t<Elements type="penta6" name="%s">.*?</Elements>\n' % d, '', t, flags=re.S)
            t, n2 = re.subn(r'\t\t<SolidDomain name="%s" mat="[^"]+" />\n' % d, '', t)
            t, n3 = re.subn(r'\t*<ElementData type="mat_axis" elem_set="%s">.*?</ElementData>\n' % d, '', t, flags=re.S)
            assert (n1, n2, n3) == (1, 1, 1), (d, n1, n2, n3)
    for c in CONTACTS:
        pr = re.search(r'\t\t<SurfacePair name="%s">\s*<primary>([^<]+)</primary>\s*<secondary>([^<]+)</secondary>\s*</SurfacePair>\n' % c, t)
        assert pr, c
        t = t[:pr.start()] + t[pr.end():]
        for s in pr.groups():
            t, n = re.subn(r'\t\t<Surface name="%s">.*?</Surface>\n' % s, '', t, flags=re.S)
            assert n == 1 and f'"{s}"' not in t and f'>{s}<' not in t, s
        t, n = re.subn(r'\t\t<contact name="%s" [^>]*>.*?</contact>\n' % c, '', t, flags=re.S)
        assert n == 1, c
    t, n = re.subn(r'\t\t<constraint name="strap_wide_ties"[^>]*>.*?</constraint>\n', '', t, flags=re.S)
    assert 'PM_avw_bottom_left_sheet_solid' not in t and 'sheet_wide' not in t
    msg = ('NOT IN SOURCE (spring layout): the fan straps replaced by single springs like PM_conn (strap edge on the curve / AVW -> the '
           f'strap\'s fixed PM-arch anchors, in order): ' + '; '.join(report) + f'; material PM_avw_strap_conn_mat = PM_conn_mat\'s table, '
           f'scale {scale:.4g} (= {total:g} PM_conn springs a side); strap elements / domains / mat_axis, contacts {", ".join(CONTACTS)} and '
           f'strap_wide_ties removed ({n} tie block). One change.')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch454.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--total') + 1]) if '--total' in a else 5 / 3)
