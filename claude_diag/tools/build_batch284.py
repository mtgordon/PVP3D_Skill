"""Batch 284 (2026-10-14; the user: L283_junction_fine2z18_MAX "the sides still penetrate themselves around 0.65 on one side and
0.69 on the other"; "go ahead with the search radius 1 mm MAX run ... focus on eliminating the penetration").
Measured (curve surface nodes inside another curve hex, none at rest): both curves fold into themselves at their distal (PeB)
ends in the split zone (rest z -28..-41); 0.36 mm deep at t 0.65, 0.52 at 0.66, 1.1 mm by 0.68 (to 1.36 at 0.8). The curves'
self-contact curve_all_self is one-pass PENALTY (penalty 5, auto) with search_radius 0.5 mm: past 0.5 mm it loses the node.
One change per variant to the curve_all_self contact block only (NOT IN SOURCE: the self-contact itself is ours):
  sr10    search_radius 1.0   (the user's request; gotcha 30: 1.0 crawled from t 0 on an older curve version)
  sr075   search_radius 0.75  (below the 0.86 mm smallest rest gap to a non-neighbour curve facet)
  pen25   penalty 25 (x5)
  pen100  penalty 100 (x20)
  auglag  laugon AUGLAG (augmented Lagrangian, tolerance 0.1, maxaug 15 as set)
usage: py -3.10 build_batch284.py SRC_RUN VARIANT [VARIANT ...] [--prefix L284_junction]  (case taken from SRC_RUN's suffix)"""
import os
import re
import sys

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
VARIANTS = {'sr10': ('search_radius', '1.0'), 'sr075': ('search_radius', '0.75'), 'pen25': ('penalty', '25'),
            'pen100': ('penalty', '100'), 'auglag': ('laugon', 'AUGLAG')}


def main(argv):
    prefix = 'L284_junction'
    if '--prefix' in argv:
        i = argv.index('--prefix'); prefix = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    src, variants = argv[0], argv[1:]
    case = src.rsplit('_', 1)[1]
    txt = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    m = re.search(r'<contact name="curve_all_self".*?</contact>', txt, re.S)
    assert m, 'no curve_all_self contact'
    for v in variants:
        tag, val = VARIANTS[v]
        blk, n = re.subn(rf'<{tag}>[^<]*</{tag}>', f'<{tag}>{val}</{tag}>', m.group(0))
        assert n == 1, (tag, n)
        new = f'{prefix}_{v}_{case}'
        d = os.path.join(RUNS, new); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, new + '.feb'), 'w', encoding='utf-8', newline='') as f:
            f.write(txt[:m.start()] + blk + txt[m.end():])
        with open(os.path.join(d, new + '.feb.changes.txt'), 'w', encoding='utf-8') as f:
            f.write(f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch284.py)\n'
                    f'curve_all_self (the curves\' self-contact, NOT IN SOURCE): {tag} -> {val} (one change)\n')
        print(new, tag, val)


if __name__ == '__main__':
    main(sys.argv[1:])
