"""Batch 310 (2026-10-17; the user saw the curves wrinkle into folds along their length (pictures of L294b) and said "try both":
(1) a pre-bend of the curve end, (2) a contact designed for the folds; this is (2)).
NOT IN SOURCE: the curve_all facets that fold in L294b MAX at t 0.74-0.78 (curve_selfpen_study/fold_zone.py: facets with >= 2
nodes among the penetrating nodes, the hexes they go into, and hexes shortened along the canal below 0.8; 452 of 1632) move to
their own self-contact curve_fold_self (a copy of curve_all_self) with offset OFF (below the smallest 0.61 mm rest gap between
non-neighbouring facets of the set), so a fold is pushed apart before it closes; the other facets keep curve_all_self.
usage: py -3.10 build_batch310.py SRC_RUN NEW_RUN NODEFILE [--off 0.4] [--pen P]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new, nodefile, off, pen):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    S = {int(v) for v in open(nodefile).read().split()}
    ms = re.search(r'(\t\t<Surface name="curve_all">\n)(.*?)(\t\t</Surface>\n)', t, re.S)
    rows = re.findall(r'\t\t\t<(\w+) id="\d+">([0-9,]+)</\1>\n', ms.group(2))
    keep, fold = [], []
    for tag, c in rows:
        (fold if sum(int(v) in S for v in c.split(',')) >= 2 else keep).append((tag, c))
    body = lambda L: ''.join(f'\t\t\t<{g} id="{k}">{c}</{g}>\n' for k, (g, c) in enumerate(L, 1))   # noqa: E731
    t = t.replace(ms.group(0), ms.group(1) + body(keep) + ms.group(3) + '\t\t<Surface name="curve_fold">\n' + body(fold) + ms.group(3), 1)
    mp = re.search(r'\t\t<SurfacePair name="curve_all_self">.*?</SurfacePair>\n', t, re.S).group(0)
    t = t.replace(mp, mp + mp.replace('curve_all_self', 'curve_fold_self').replace('>curve_all<', '>curve_fold<'), 1)
    mc = re.search(r'\t\t<contact name="curve_all_self".*?</contact>\n', t, re.S).group(0)
    c2 = mc.replace('name="curve_all_self" surface_pair="curve_all_self"', 'name="curve_fold_self" surface_pair="curve_fold_self"')
    c2 = re.sub(r'<offset>[^<]*</offset>', f'<offset>{off:g}</offset>', c2) if '<offset>' in c2 else c2.replace('</contact>', f'\t<offset>{off:g}</offset>\n\t\t</contact>')
    if pen is not None:
        c2 = re.sub(r'<penalty>[^<]*</penalty>', f'<penalty>{pen:g}</penalty>', c2)
    t = t.replace(mc, mc + c2, 1)
    assert 'curve_fold_self' in t and f'<offset>{off:g}</offset>' in c2
    msg = (f'NOT IN SOURCE: {len(fold)} of {len(rows)} curve_all facets (the folding zone from {os.path.basename(nodefile)}) moved to '
           f'curve_fold with its own self-contact curve_fold_self (copy of curve_all_self, offset {off:g}' +
           (f', penalty {pen:g}' if pen is not None else '') + ')')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch310.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], a[2], g('--off', 0.4), g('--pen', None))
