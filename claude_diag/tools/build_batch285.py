"""Batch 285 (2026-10-14; follows batch 284). Batch 284 changed the curves' whole self-contact (curve_all_self) to stop the
distal fold-through in L283_junction_fine2z18_MAX: every variant (search radius 0.75 / 1.0, penalty x5 / x20, AUGLAG) stuck
at t 0.0003-0.008 (control L284_junction_ctrl_MAX = L283 ran normally). The stuck spot is NOT the distal fold: it is the
curves' upper end, where canal_wrap_hex_avw_cervix meets apex_contour (rest about (+-20, -30, +31); curve_all holds 24 + 24
apex_contour faces, and 2 faces there are in contact from t 0.005 in L283). A wider / stiffer self-contact there pushes the
two meshes apart (to 15 mm by t 0.004).
Here: curve_all's facets with every node at rest z < ZCUT (-15 mm; the fold is at rest z -28..-41) move to a new surface
curve_distal with its own self-contact curve_distal_self (a copy of curve_all_self); curve_all keeps the rest, with its contact
unchanged. Then ONE change in curve_distal_self per variant (NOT IN SOURCE: both self-contacts are ours):
  split   no change (the split alone; contact between the two zones across z -15 is lost)
  sr06    search_radius 0.6
  sr075   search_radius 0.75 (the smallest rest gap to a non-neighbour curve facet there is 0.86 mm)
  sr10    search_radius 1.0 (the user's request; 24 distal nodes are 0.86-1.0 mm from such a facet at rest)
  pen7    penalty 7.5 (x1.5)
  pen10   penalty 10 (x2)
  pen15   penalty 15 (x3)
  pen25   penalty 25 (x5)
  pen50   penalty 50 (x10)
  pen100  penalty 100 (x20)
  auglag  laugon AUGLAG
  tp      two_pass 1
  off03 / off06  offset 0.3 / 0.6 mm (the faces are held that far apart; tested on 2 blocks 10-15: radius need not grow;
          an offset larger than a rest gap errors at t 0 -> must stay < 0.86 mm, the smallest distal rest gap)
  seg0    seg_up 0 (contact pairing updated every iteration; 2 = frozen after two updates per step)
  snf     contact type sliding-node-on-facet (closest-point projection; sliding-elastic casts along the node normal, which
          can miss the opposite face in a closing crease)
usage: py -3.10 build_batch285.py SRC_RUN VARIANT [VARIANT ...] [--zcut -15] [--prefix L285_distal]"""
import os
import re
import sys

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
VARIANTS = {'split': None, 'sr06': ('search_radius', '0.6'), 'sr075': ('search_radius', '0.75'), 'sr10': ('search_radius', '1.0'),
            'pen7': ('penalty', '7.5'), 'pen10': ('penalty', '10'), 'pen15': ('penalty', '15'), 'pen25': ('penalty', '25'), 'pen50': ('penalty', '50'), 'pen100': ('penalty', '100'), 'auglag': ('laugon', 'AUGLAG'), 'tp': ('two_pass', '1'), 'off03': ('offset', '0.3'), 'off06': ('offset', '0.6'), 'off045': ('offset', '0.45'), 'off06tp': [('offset', '0.6'), ('two_pass', '1')], 'seg0': ('seg_up', '0'), 'snf': ('type', 'sliding-node-on-facet')}


def surface_block(name, facets):
    rows = ''.join(f'\t\t\t<{t} id="{i}">{c}</{t}>\n' for i, (t, c) in enumerate(facets, 1))
    return f'\t\t<Surface name="{name}">\n{rows}\t\t</Surface>\n'


def main(argv):
    prefix, zcut = 'L285_distal', -15.0
    for opt in ('--prefix', '--zcut'):
        if opt in argv:
            i = argv.index(opt); v = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
            if opt == '--prefix':
                prefix = v
            else:
                zcut = float(v)
    src, variants = argv[0], argv[1:]
    case = src.rsplit('_', 1)[1]
    txt = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    z = {int(m.group(1)): float(m.group(2).split(',')[2])
         for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', txt)}
    ms = re.search(r'\t\t<Surface name="curve_all">\n(.*?)\t\t</Surface>\n', txt, re.S)
    facets = re.findall(r'<(\w+) id="\d+">([^<]+)</\1>', ms.group(1))
    prox, dist = [], []
    for t, c in facets:
        (dist if all(z[int(n)] < zcut for n in c.split(',')) else prox).append((t, c))
    print(f'curve_all {len(facets)} facets -> curve_all {len(prox)} + curve_distal {len(dist)} (all nodes at rest z < {zcut})')
    mp = re.search(r'(\t\t<SurfacePair name="curve_all_self">.*?</SurfacePair>\n)', txt, re.S)
    mc = re.search(r'(\t\t<contact name="curve_all_self".*?</contact>\n)', txt, re.S)
    assert ms and mp and mc
    pair = mp.group(1).replace('curve_all_self', 'curve_distal_self').replace('curve_all', 'curve_distal')
    for v in variants:
        blk = mc.group(1).replace('curve_all_self', 'curve_distal_self')
        if v == 'snf':
            # sliding-node-on-facet keeps every sliding-elastic value it accepts (smooth_aug, symmetric_stiffness, tension,
            # flip_*, shell_bottom_*, offset dropped; smooth_aug tested: rejected)
            keep = ('laugon', 'tolerance', 'gaptol', 'penalty', 'auto_penalty', 'update_penalty', 'two_pass', 'knmult',
                    'search_tol', 'search_radius', 'seg_up', 'minaug', 'maxaug', 'node_reloc', 'fric_coeff')
            rows = ''.join(m.group(0) + '\n' for m in re.finditer(r'\t\t\t<(\w+)>[^<]*</\1>', blk) if m.group(1) in keep)
            blk = ('\t\t<contact name="curve_distal_self" surface_pair="curve_distal_self" type="sliding-node-on-facet">\n'
                   + rows + '\t\t</contact>\n')
        elif VARIANTS[v]:
            for tag, val in (VARIANTS[v] if isinstance(VARIANTS[v], list) else [VARIANTS[v]]):
                blk, n = re.subn(rf'<{tag}>[^<]*</{tag}>', f'<{tag}>{val}</{tag}>', blk)
                assert n == 1
        out = txt[:ms.start()] + surface_block('curve_all', prox) + surface_block('curve_distal', dist) + txt[ms.end():]
        out = out.replace(mp.group(1), mp.group(1) + pair, 1)
        out = out.replace(mc.group(1), mc.group(1) + blk, 1)
        new = f'{prefix}_{v}_{case}'
        d = os.path.join(RUNS, new); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, new + '.feb'), 'w', encoding='utf-8', newline='') as f:
            f.write(out)
        chs = VARIANTS[v] if isinstance(VARIANTS[v], list) else [VARIANTS[v]] if VARIANTS[v] else []
        ch = ', '.join(f'{a} -> {b}' for a, b in chs) or 'none (the split alone)'
        with open(os.path.join(d, new + '.feb.changes.txt'), 'w', encoding='utf-8') as f:
            f.write(f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch285.py)\n'
                    f'NOT IN SOURCE: curve_all facets with every node at rest z < {zcut} ({len(dist)} of {len(facets)}) moved '
                    f'to curve_distal with its own self-contact curve_distal_self (copy of curve_all_self); '
                    f'change in curve_distal_self: {ch}\n')
        print(new, ch)


if __name__ == '__main__':
    main(sys.argv[1:])
