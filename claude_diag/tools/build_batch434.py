"""Batch 434 (2026-10-23; the user: "the vaginal wall would not be stiff to compression as it would fold"; strip test
shell_tests/curve_fibre/gen.py -> "probably run the fiber one next"). NOT IN SOURCE (the curves are ours): the curves' material
Curve_soft (Yeoh c1 0.00125, c2 0.00396, k 0.21; in a 1-hex uniaxial test 3.00 kPa at +30 % and 6.94 kPa at -30 %: stiffer in
compression) is replaced in canal_wrap_hex_avw / _pvw / _pvw_fine by Curve_fibre: an uncoupled solid mixture of a soft
Mooney-Rivlin matrix (c1 C1M = 0.1 x Curve_soft c1, c2 0), k K, and two tension-only fiber-pow-linear-uncoupled fibres (E, beta,
lam0) along the local material axes: a = along the canal (the hex's axial edges 0-4, 1-5, 2-6, 3-7 averaged), d = around the U
(edges 0-1, 3-2, 4-5, 7-6 averaged). Strip test: tension 0.42 / 0.96 / 1.49 / 2.02 / 2.54 kPa at +6..+30 % (Curve_soft 0.44 /
0.91 / 1.46 / 2.15 / 3.00), compression 0.50 kPa at -30 % (14x softer). The curve parts on the cervix (canal_wrap_hex_avw_cervix,
Vagina_Cervix) are not changed. Mass damping and density unchanged.
usage: py -3.10 build_batch434.py SRC_RUN NEW_RUN [--c1m 0.000125] [--e 0.0085] [--beta 2.5] [--lam0 1.04] [--k 0.21]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

DOMS = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')


def main(src, new, c1m, E, beta, lam0, k):
    t, X, E_ = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    fib = lambda v: (f'\t\t\t<solid type="fiber-pow-linear-uncoupled">\n\t\t\t\t<fiber type="vector">{v}</fiber>\n'   # noqa: E731
                     f'\t\t\t\t<E>{E:g}</E>\n\t\t\t\t<beta>{beta:g}</beta>\n\t\t\t\t<lam0>{lam0:g}</lam0>\n\t\t\t</solid>\n')
    mat = (f'\t\t<material id="{mid}" name="Curve_fibre" type="uncoupled solid mixture">\n\t\t\t<density>1.06e-09</density>\n'
           f'\t\t\t<k>{k:g}</k>\n\t\t\t<solid type="Mooney-Rivlin">\n\t\t\t\t<c1>{c1m:g}</c1>\n\t\t\t\t<c2>0</c2>\n\t\t\t</solid>\n'
           + fib('1,0,0') + fib('0,1,0') + '\t\t</material>\n')
    i = t.index('</Material>')
    t = t[:i] + mat + '\t' + t[i:]
    data, worst = '', 1.0
    doms = re.findall(r'<SolidDomain name="(canal_wrap_hex[^"]*)" mat="Curve_soft" />', t)   # 10-23: + split domains (_len, _around)
    assert set(DOMS) <= set(doms), doms
    for d in doms:
        t = t.replace(f'<SolidDomain name="{d}" mat="Curve_soft" />', f'<SolidDomain name="{d}" mat="Curve_fibre" />')
        els = [c for n, c in E_.values() if n == d]
        data += f'\t\t<ElementData type="mat_axis" elem_set="{d}">\n'
        for lid, c in enumerate(els, 1):
            P = np.array([X[v] for v in c])
            a = sum(P[q + 4] - P[q] for q in range(4)); a /= np.linalg.norm(a)
            u = (P[1] - P[0]) + (P[2] - P[3]) + (P[5] - P[4]) + (P[6] - P[7])
            u = u - (u @ a) * a; u /= np.linalg.norm(u)
            worst = min(worst, np.linalg.norm(u))
            data += f'\t\t\t<elem lid="{lid}"><a>{a[0]:.6f},{a[1]:.6f},{a[2]:.6f}</a><d>{u[0]:.6f},{u[1]:.6f},{u[2]:.6f}</d></elem>\n'
        data += '\t\t</ElementData>\n'
    if '</MeshData>' in t:
        i = t.index('</MeshData>'); t = t[:i] + data + '\t' + t[i:]
    else:
        i = t.index('</MeshDomains>') + len('</MeshDomains>\n'); t = t[:i] + '\t<MeshData>\n' + data + '\t</MeshData>\n' + t[i:]
    msg = (f'NOT IN SOURCE (our curves): canal_wrap_hex_avw / _pvw / _pvw_fine Curve_soft (Yeoh) -> Curve_fibre (uncoupled mixture: '
           f'Mooney-Rivlin c1 {c1m:g}, k {k:g}, + 2 tension-only fiber-pow-linear-uncoupled E {E:g} beta {beta:g} lam0 {lam0:g} along '
           f'the canal (a) and around the U (d), per-element mat_axis); strip test: tension within 5-24 % of Curve_soft to +30 %, '
           f'compression 14x softer')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch434.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: float(a[a.index(k_) + 1]) if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--c1m', 0.000125), g('--e', 0.0085), g('--beta', 2.5), g('--lam0', 1.04), g('--k', 0.21))
