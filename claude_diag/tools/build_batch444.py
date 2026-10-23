"""Batch 444 (2026-10-23; L442c review). NOT IN SOURCE (wall material): L442c's Wall_fibre (2 in-plane fibres at 0 / 90 deg of
the element's a / d axes, k 0.035) let the AVW stretch diagonally and change volume (V/V0 0.48-1.59 at t 0.25, exploded at
t 0.317; L440g with the Yeoh walls stays 0.90-1.09). Strip test (shell_tests/curve_fibre/gen.py mixa): pulled at 45 deg
the 2-fibre wall is 6x softer than along a fibre (0.79 vs 3.91 kPa at +30 %; source Yeoh 5.01). The user (10-23): the tissue
should be soft in compression but stiff in tension. This rewrites Wall_fibre as 4 tension-only fibres at 0 / 45 / 90 / 135 deg
in the wall plane (the same per-element mat_axis), E 0.014, beta 2.5, lam0 1.04, matrix c1 0.000208, and the source wall's
k 0.35 (k does not stiffen uniaxial compression; it stops the volume drift). Strip test, E 0.014 k 0.35: tension at 0 deg
0.71 / 1.64 / 2.60 / 3.61 / 4.69 kPa at +6..+30 %, at 22.5 deg 0.65 / 1.71 / 2.88 / 4.12 / 5.41 (Yeoh 0.73 / 1.51 / 2.44 /
3.59 / 5.01); +60 %: 11.3 (Yeoh 17.0); compression -30 %: 0.83 (Yeoh 11.57, 14x softer).
usage: py -3.10 build_batch444.py SRC_RUN NEW_RUN [--e 0.014] [--k 0.35] [--angles 0,45,90,135] [--c1m 0.000208]
L447 (10-23): --c1m 0.0005 --e 0.011 (L444/L445 AVW V/V0 still 0.77-1.24 at t 0.215 vs Yeoh 0.94-1.06): tension at +30 %
4.37 / 4.84 kPa (0 / 22.5 deg), compression -30 % 1.86 kPa (6x softer than Yeoh)."""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, E, k, angs, c1m=0.000208, beta=2.5, lam0=1.04):
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8').read()
    m = re.search(r'(\t*)<material id="(\d+)" name="Wall_fibre" type="uncoupled solid mixture">.*?</material>\n', t, re.S)
    assert m, 'no Wall_fibre in ' + src
    ind = m.group(1)
    fib = ''.join(f'{ind}\t<solid type="fiber-pow-linear-uncoupled">\n{ind}\t\t<fiber type="vector">{math.cos(math.radians(a)):.6f},'
                  f'{math.sin(math.radians(a)):.6f},0</fiber>\n{ind}\t\t<E>{E:g}</E>\n{ind}\t\t<beta>{beta:g}</beta>\n'
                  f'{ind}\t\t<lam0>{lam0:g}</lam0>\n{ind}\t</solid>\n' for a in angs)
    mat = (f'{ind}<material id="{m.group(2)}" name="Wall_fibre" type="uncoupled solid mixture">\n{ind}\t<density>1.06e-09</density>\n'
           f'{ind}\t<k>{k:g}</k>\n{ind}\t<solid type="Mooney-Rivlin">\n{ind}\t\t<c1>{c1m:g}</c1>\n{ind}\t\t<c2>0</c2>\n{ind}\t</solid>\n'
           + fib + f'{ind}</material>\n')
    t = t[:m.start()] + mat + t[m.end():]
    msg = (f'NOT IN SOURCE (wall material): Wall_fibre (AVW / PVW) -> {len(angs)} in-plane tension-only fiber-pow-linear-uncoupled '
           f'fibres at {",".join(f"{a:g}" for a in angs)} deg of the per-element mat_axis, E {E:g} beta {beta:g} lam0 {lam0:g}, matrix '
           f'Mooney-Rivlin c1 {c1m:g}, k {k:g} (source wall k); strip test: tension within ~10 % of the wall law to +30 % in any '
           f'in-plane direction, compression 14x softer (E 0.014, c1m 0.000208; see the docstring for other values)')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch444.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    opt = lambda n, d: (a[a.index(n) + 1] if n in a else d)  # noqa: E731
    main(a[0], a[1], float(opt('--e', 0.014)), float(opt('--k', 0.35)), [float(v) for v in opt('--angles', '0,45,90,135').split(',')],
         float(opt('--c1m', 0.000208)))
