"""Batch 442 (2026-10-23; the user: "include the avw and pvw in the fiber modification as well. Let's do that for future
versions"). NOT IN SOURCE (the source walls are Yeoh-fitted Marlow, Vagina_AVW / Vagina_PVW c1 0.00209, c2 0.00659, k 0.35;
1-hex strip test: 5.01 kPa at +30 %, 11.57 kPa at -30 %: stiffer in compression). The AVW (_PickedSet347) and PVW
(_PickedSet64) become Wall_fibre (one tissue, one material, as the source's two identical laws): an uncoupled solid mixture
of a soft Mooney-Rivlin matrix (c1 C1M) + two tension-only fiber-pow-linear-uncoupled fibres (E, beta, lam0) in the wall's
plane, k K (scaled with the matrix as the curves' L436). In-plane axes per element: every AVW / PVW hex has its through-
thickness direction along its 0-4 edges (checked 10-23 against the nearest lumen face: 3796 of 3796, median margin 0.87-0.92),
so a = the 0-1 edge family, d = the 0-3 family orthogonalised against a. Strip test (shell_tests/curve_fibre/gen.py
mix 0.000208 0.01417 2.5 1.04 0.035): tension 0.67 / 1.51 / 2.33 / 3.13 / 3.91 kPa at +6..+30 % (Yeoh 0.73 / 1.51 / 2.44 /
3.59 / 5.01), compression 0.83 kPa at -30 % (14x softer). The curves are converted separately (build_batch434).
usage: py -3.10 build_batch442.py SRC_RUN NEW_RUN [--c1m 0.000208] [--e 0.01417] [--beta 2.5] [--lam0 1.04] [--k 0.035]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

WALLS = {'_PickedSet347': 'Vagina_AVW', '_PickedSet64': 'Vagina_PVW'}


def main(src, new, c1m, E, beta, lam0, k):
    t, X, E_ = load(src)
    X = {n: np.asarray(p, float) for n, p in dict(X).items()}
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    fib = lambda v: (f'\t\t\t<solid type="fiber-pow-linear-uncoupled">\n\t\t\t\t<fiber type="vector">{v}</fiber>\n'   # noqa: E731
                     f'\t\t\t\t<E>{E:g}</E>\n\t\t\t\t<beta>{beta:g}</beta>\n\t\t\t\t<lam0>{lam0:g}</lam0>\n\t\t\t</solid>\n')
    mat = (f'\t\t<material id="{mid}" name="Wall_fibre" type="uncoupled solid mixture">\n\t\t\t<density>1.06e-09</density>\n'
           f'\t\t\t<k>{k:g}</k>\n\t\t\t<solid type="Mooney-Rivlin">\n\t\t\t\t<c1>{c1m:g}</c1>\n\t\t\t\t<c2>0</c2>\n\t\t\t</solid>\n'
           + fib('1,0,0') + fib('0,1,0') + '\t\t</material>\n')
    i = t.index('</Material>')
    t = t[:i] + mat + '\t' + t[i:]
    data = ''
    for d, old in WALLS.items():
        tag = f'<SolidDomain name="{d}" mat="{old}" />'
        assert tag in t, tag
        t = t.replace(tag, f'<SolidDomain name="{d}" mat="Wall_fibre" />')
        els = [c for n, c in E_.values() if n == d and len(c) == 8]
        data += f'\t\t<ElementData type="mat_axis" elem_set="{d}">\n'
        for lid, c in enumerate(els, 1):
            P = np.array([X[v] for v in c])
            a = (P[1] - P[0]) + (P[2] - P[3]) + (P[5] - P[4]) + (P[6] - P[7]); a /= np.linalg.norm(a)
            u = (P[3] - P[0]) + (P[2] - P[1]) + (P[7] - P[4]) + (P[6] - P[5])
            u = u - (u @ a) * a; u /= np.linalg.norm(u)
            data += f'\t\t\t<elem lid="{lid}"><a>{a[0]:.6f},{a[1]:.6f},{a[2]:.6f}</a><d>{u[0]:.6f},{u[1]:.6f},{u[2]:.6f}</d></elem>\n'
        data += '\t\t</ElementData>\n'
    if '</MeshData>' in t:
        i = t.index('</MeshData>'); t = t[:i] + data + '\t' + t[i:]
    else:
        i = t.index('</MeshDomains>') + len('</MeshDomains>\n'); t = t[:i] + '\t<MeshData>\n' + data + '\t</MeshData>\n' + t[i:]
    msg = (f'NOT IN SOURCE (wall material): AVW / PVW Vagina_AVW / Vagina_PVW (Yeoh) -> Wall_fibre (uncoupled mixture: Mooney-Rivlin '
           f'c1 {c1m:g}, k {k:g}, + 2 in-plane tension-only fiber-pow-linear-uncoupled E {E:g} beta {beta:g} lam0 {lam0:g}, per-element '
           f'mat_axis, thickness along the hexes\' 0-4 edges); strip test: tension within 8-22 % of the wall law to +30 %, compression 14x softer')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch442.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: float(a[a.index(k_) + 1]) if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--c1m', 0.000208), g('--e', 0.01417), g('--beta', 2.5), g('--lam0', 1.04), g('--k', 0.035))
