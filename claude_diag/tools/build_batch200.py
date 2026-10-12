"""Batch 200 (2026-10-07, the user: "yes, start the PVW-side softening tube test"). On the tube the AVW and PVW descend
together (joined all round through the curves): posterior x0.4 moved Bp +9.5 but Ba +6.9 (L197). To let the PVW side give
relative to the AVW, one LOCAL change on L192_botfull_distpress_cdamp1000_P1 (B: equal PM-arch springs, distal-half curve
pressure, curve damping C 1000; the smallest fold), P1:
  L200_distpress_pvwcurve05_P1   the PVW half of the curves (domain canal_wrap_hex_pvw, 1320 elements) gets its own material
                                 Vagina_PVW_curve = Vagina_PVW c1, c2, k x 0.5 (NOT IN SOURCE; the PVW wall itself unchanged)
usage: py -3.10 build_batch200.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L192_botfull_distpress_cdamp1000_P1'
NAME = 'L200_distpress_pvwcurve05_P1'
FACTOR = 0.5

if __name__ == '__main__':
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    a = t.rindex('<material', 0, t.index('name="Vagina_PVW" type="Yeoh">'))
    b = t.index('</material>', a) + len('</material>')
    blk = t[a:b]
    for k in ('c1', 'c2', 'k'):
        blk = re.sub(rf'<{k}>([^<]+)</{k}>', lambda m: f'<{k}>{float(m.group(1)) * FACTOR:.7g}</{k}>', blk)
    nid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    blk = re.sub(r'id="\d+"', f'id="{nid}"', blk, count=1).replace('name="Vagina_PVW"', 'name="Vagina_PVW_curve"')
    i = t.index('\t</Material>')
    t = t[:i] + '\t\t' + blk + '\n' + t[i:]
    old = '<SolidDomain name="canal_wrap_hex_pvw" mat="Vagina_PVW" />'
    assert t.count(old) == 1
    t = t.replace(old, '<SolidDomain name="canal_wrap_hex_pvw" mat="Vagina_PVW_curve" />')
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch200.py)\n'
                   f'NOT IN SOURCE (LOCAL: the PVW half of the tube curves only): canal_wrap_hex_pvw -> material Vagina_PVW_curve = '
                   f'Vagina_PVW c1, c2, k x {FACTOR:g} (walls x0.35 already, so x{0.35 * FACTOR:g} of the source); the PVW wall unchanged\n')
