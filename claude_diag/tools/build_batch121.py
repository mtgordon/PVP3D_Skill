"""Batch 121 (2026-09-29, the user's request: get the tube_wrap to open more): L119_tube_wrap with a softer wrap.
One change each from L119_tube_wrap (NOT IN SOURCE): the wrap's two domains (canal_wrap_hex / canal_wrap_wedge) get their
own copy of Vagina_AVW (Yeoh) with c1, c2 and k all scaled by F (same Poisson ratio), so the fold bends open more easily.
Tissue joins ~0.035-0.31 N/mm per seam node vs the 0.01 N/mm seam springs, so F ~0.03-0.1 puts the wrap near the springs.
  L121_tube_wrap_soft030 / 010 / 003 / 001 / 0003   F = 0.3 / 0.1 / 0.03 / 0.01 / 0.003
usage: py -3.10 build_batch121.py"""
import os
import re

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
SRC = 'L119_tube_wrap'
C1, C2, K = 0.005964127, 0.01883547, 1.0
VARIANTS = {'soft030': 0.3, 'soft010': 0.1, 'soft003': 0.03, 'soft001': 0.01, 'soft0003': 0.003}

txt = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
mid = max(int(x) for x in re.findall(r'<material id="(\d+)"', txt)) + 1
for tag, F in VARIANTS.items():
    name = f'L121_tube_wrap_{tag}'
    d = os.path.join(RUNS, name)
    assert not os.path.exists(d), f'{name} exists; not overwriting'
    mat = (f'\t\t<material id="{mid}" name="Canal_wrap_soft" type="Yeoh">\n\t\t\t<density>1.06e-09</density>\n'
           f'\t\t\t<c1>{C1 * F:.6g}</c1>\n\t\t\t<c2>{C2 * F:.6g}</c2>\n\t\t\t<k>{K * F:.6g}</k>\n\t\t</material>\n')
    t = txt.replace('\t</Material>', mat + '\t</Material>', 1)
    assert t != txt
    n0 = t.count('mat="Vagina_AVW" />')
    t = t.replace('<SolidDomain name="canal_wrap_hex" mat="Vagina_AVW" />', '<SolidDomain name="canal_wrap_hex" mat="Canal_wrap_soft" />')
    t = t.replace('<SolidDomain name="canal_wrap_wedge" mat="Vagina_AVW" />', '<SolidDomain name="canal_wrap_wedge" mat="Canal_wrap_soft" />')
    assert t.count('mat="Canal_wrap_soft"') == 2, n0
    os.makedirs(d)
    open(os.path.join(d, name + '.feb'), 'w', encoding='utf-8', newline='').write(t)
    open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8').write(
        f'{name}: {SRC} + the wrap (canal_wrap_hex / canal_wrap_wedge) in its own Yeoh, Vagina_AVW x {F} '
        f'(c1 {C1 * F:.6g}, c2 {C2 * F:.6g}, k {K * F:.6g}); NOT IN SOURCE\n')
    print(name, 'written, F =', F)
