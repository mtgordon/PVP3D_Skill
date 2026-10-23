"""Batch 420 (2026-10-22; L413b: the curve pokes through the strap; 2-6 crossings a step through the strap's BOTTOM / SIDE faces,
which are not in curve_strap_secondary (top faces only)). NOT IN SOURCE (our contact): every exterior face of the strap domains
(PM_avw_bottom_*_sheet_solid / _sheet_wide) is put in curve_strap_secondary, outward (bottom (0,2,1), top (3,4,5), sides
(a,b,b+3,a+3)), except faces with a node shared with the curves (the tied rim). Faces already there are kept as they are.
usage: py -3.10 build_batch420.py SRC_RUN NEW_RUN [--rim]
--rim (10-23, L435: the first curve nodes through the strap went through the side / bottom faces AT the tied rim): include
the faces touching curve nodes too (no curve node within 0.6 mm of them at rest)"""
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402

PF = [((0, 2, 1), 'bottom'), ((3, 4, 5), 'top'), ((0, 1, 4, 3), 'side'), ((1, 2, 5, 4), 'side'), ((2, 0, 3, 5), 'side')]


def main(src, new, rim=False):
    t, X, E = load(src)
    W = [c for n, c in E.values() if n.startswith('PM_avw_bottom_') and len(c) == 6]
    curve = {v for n, c in E.values() if n.startswith('canal_wrap') for v in c}
    cnt = Counter(frozenset(c[i] for i in f) for c in W for f, _ in PF)
    m = re.search(r'(<Surface name="curve_strap_secondary">\n)(.*?)(\t\t</Surface>)', t, re.S)
    have = [(k, list(map(int, f.split(',')))) for k, f in re.findall(r'<(tri3|quad4) id="\d+">([\d,]+)</', m.group(2))]
    hs = {frozenset(f) for _, f in have}
    add, ty = [], Counter()
    for c in W:
        for f, name in PF:
            nodes = [c[i] for i in f]
            s = frozenset(nodes)
            if cnt[s] != 1 or s in hs or (s & curve and not rim):
                continue
            add.append(('tri3' if len(nodes) == 3 else 'quad4', nodes)); ty[name] += 1; hs.add(s)
    body = ''.join(f'\t\t\t<{k} id="{i}">' + ','.join(map(str, f)) + f'</{k}>\n'
                   for i, (k, f) in enumerate(have + add, 1))
    t = t[:m.start(2)] + body + t[m.end(2):]
    msg = (f'NOT IN SOURCE (contact surface): curve_strap_secondary {len(have)} -> {len(have) + len(add)} facets: + {dict(ty)} '
           f'(all exterior strap faces, outward; faces touching curve nodes ' + ('INCLUDED' if rim else 'left out') + ')')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch420.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], '--rim' in sys.argv)
