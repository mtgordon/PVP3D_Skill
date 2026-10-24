"""Batch 449 (2026-10-23; the user: "it seems like the fiber directions aren't working very well ... Is there a way to do tensile
stiffness without much compressive stiffness without doing fibers?"). NOT IN SOURCE (wall material): the AVW / PVW Wall_fibre
(4 in-plane tension-only fibres; L444 / L445: AVW distal-corner elements drift and fold) -> Wall_ogden, an isotropic 2-term
uncoupled Ogden fitted to the source wall's uniaxial TENSION (Yeoh c1 0.00209, c2 0.00659): c1 0.00221 m1 2, c2 0.00563 m2 8,
k 0.35 (source k). FEBio strip test (shell_tests/curve_fibre/gen.py ogden): tension 0.72 / 1.50 / 2.42 / 3.54 / 4.93 kPa at
+6..+30 % (Yeoh 0.73 / 1.51 / 2.44 / 3.59 / 5.01), compression -30 % 5.58 kPa (Yeoh 11.57: 2.1x softer). Any isotropic law has
the same small-strain stiffness in tension and compression; the asymmetry builds up beyond ~10 % strain. The walls' mat_axis
data is removed (no fibres).
usage: py -3.10 build_batch449.py SRC_RUN NEW_RUN [--c 0.00221,0.00563] [--m 2,8] [--k 0.35] [--from Wall_fibre --to Wall_ogden]
L450 (10-23, the user): the curves too: --from Curve_fibre --to Curve_ogden --c 0.001326,0.003378 --k 0.21 (the wall fit x0.6,
as Curve_soft = 0.6 x the wall Yeoh; every domain on the material converted, its mat_axis removed)."""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, cs, ms, k, old='Wall_fibre', name='Wall_ogden'):
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8').read()
    m = re.search(r'(\t*)<material id="(\d+)" name="%s" type="uncoupled solid mixture">.*?</material>\n' % old, t, re.S)
    assert m, f'no {old} in ' + src
    doms = re.findall(r'<SolidDomain name="([^"]+)" mat="%s" />' % old, t)
    assert doms, old
    ind = m.group(1)
    terms = ''.join(f'{ind}\t<c{i}>{c:g}</c{i}>\n{ind}\t<m{i}>{mm:g}</m{i}>\n' for i, (c, mm) in enumerate(zip(cs, ms), 1))
    mat = (f'{ind}<material id="{m.group(2)}" name="{name}" type="Ogden">\n{ind}\t<density>1.06e-09</density>\n'
           f'{ind}\t<k>{k:g}</k>\n' + terms + f'{ind}</material>\n')
    t = t[:m.start()] + mat + t[m.end():]
    for d in doms:
        t = t.replace(f'<SolidDomain name="{d}" mat="{old}" />', f'<SolidDomain name="{d}" mat="{name}" />')
        t, n = re.subn(r'\t*<ElementData type="mat_axis" elem_set="%s">.*?</ElementData>\n' % d, '', t, flags=re.S)
        assert n == 1, (d, n)
    assert old not in t
    msg = (f'NOT IN SOURCE (material): {old} ({", ".join(doms)}) -> {name} (isotropic uncoupled Ogden, '
           + ', '.join(f'c{i} {c:g} m{i} {mm:g}' for i, (c, mm) in enumerate(zip(cs, ms), 1)) + f', k {k:g}; no fibres, mat_axis removed); '
           f'strip test (the wall fit; x0.6 for the curves): tension within 2 % of the Yeoh law to +30 %, compression -30 % 2.1x softer. One change.')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch449.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    opt = lambda n, d: (a[a.index(n) + 1] if n in a else d)  # noqa: E731
    main(a[0], a[1], [float(v) for v in opt('--c', '0.00221,0.00563').split(',')], [float(v) for v in opt('--m', '2,8').split(',')],
         float(opt('--k', 0.35)), opt('--from', 'Wall_fibre'), opt('--to', 'Wall_ogden'))
