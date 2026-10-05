"""Batch 123 (2026-09-29, following L122's result: the rounded lumen end opens the seam more than the shared-node wrap).
  L123_tube_round_r30            L122_tube_round_r10 / r20 with RIN = 3.0 mm (a wider U; one change: the radius)
  L123_tube_round_r20_soft010 / soft003 / soft001
                                 L122_tube_round_r20 + the wrap (canal_wrap_hex) in its own Yeoh, Vagina_AVW x 0.1 / 0.03 /
                                 0.01 (c1, c2, k; as batch 121); one change from L122_tube_round_r20
All NOT IN SOURCE. usage: py -3.10 build_batch123.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_batch122 as b122  # noqa: E402
from variants7 import Model7, emit  # noqa: E402

RUNS = b122.RUNS
C1, C2, K = 0.005964127, 0.01883547, 1.0

name = 'L123_tube_round_r30'
assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
b122.RIN = 3.0
m = Model7(os.path.join(RUNS, b122.BASE, b122.BASE + '.feb'))
_, jh, _ = b122.wrap(m)
assert (jh > 0).all()
print(name, 'hex jac min', jh.min())
emit(name, m)

SRC = 'L122_tube_round_r20'
txt = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
mid = max(int(x) for x in re.findall(r'<material id="(\d+)"', txt)) + 1
for tag, F in (('soft010', 0.1), ('soft003', 0.03), ('soft001', 0.01)):
    name = f'L123_tube_round_r20_{tag}'
    d = os.path.join(RUNS, name)
    assert not os.path.exists(d), f'{name} exists; not overwriting'
    mat = (f'\t\t<material id="{mid}" name="Canal_wrap_soft" type="Yeoh">\n\t\t\t<density>1.06e-09</density>\n'
           f'\t\t\t<c1>{C1 * F:.6g}</c1>\n\t\t\t<c2>{C2 * F:.6g}</c2>\n\t\t\t<k>{K * F:.6g}</k>\n\t\t</material>\n')
    t = txt.replace('\t</Material>', mat + '\t</Material>', 1)
    a = '<SolidDomain name="canal_wrap_hex" mat="Vagina_AVW" />'
    assert t.count(a) == 1
    t = t.replace(a, '<SolidDomain name="canal_wrap_hex" mat="Canal_wrap_soft" />')
    os.makedirs(d)
    open(os.path.join(d, name + '.feb'), 'w', encoding='utf-8', newline='').write(t)
    open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8').write(
        f'{name}: {SRC} + the wrap (canal_wrap_hex) in its own Yeoh, Vagina_AVW x {F} (c1 {C1 * F:.6g}, c2 {C2 * F:.6g}, '
        f'k {K * F:.6g}); NOT IN SOURCE\n')
    print(name, 'written')
