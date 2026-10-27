"""Batch C018 (2026-10-26, Stephen's long session). The remaining curve self-crossing after the curve band
(C010 / C011: ~0.7-0.85 mm at mid-canal, rest z -9..+1) is the curve skin passing through the curve's own END faces, the
faces whose nodes are tied into the AVW / PVW by the canal_wrap_ends linear constraints. Those faces are not in the curve
self-contact (curve_all), by design (curve_cross_diag.py). Here a separate contact catches exactly that: primary = the curve
skin (curve_all, as the self-contact), secondary = the tied end faces (new surface curve_endfaces, outward-wound from the hex
face order), sliding-elastic with the self-contact's settings (penalty 5, search radius 0.5, one pass). The ties and every
other contact are unchanged. NOT IN SOURCE.
usage: py -3.10 build_batchC018.py SRC_RUN NEW_RUN"""
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from paths import RUNS_DIR  # noqa: E402

FI = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def convert(t):
    H = []
    for m in re.finditer(r'<Elements type="hex8" name="(canal_wrap_hex_(?!avw_cervix)[^"]*)">(.*?)</Elements>', t, re.S):
        H += [[int(v) for v in e.split(',')] for e in re.findall(r'<elem id="\d+">([^<]+)</elem>', m.group(2))]
    cnt = Counter(tuple(sorted(h[i] for i in f)) for h in H for f in FI)
    ext = [[h[i] for i in f] for h in H for f in FI if cnt[tuple(sorted(h[i] for i in f))] == 1]
    cblk = re.search(r'<constraint name="canal_wrap_ends" type="linear constraint">(.*?)</constraint>', t, re.S).group(1)
    tied = {int(v) for v in re.findall(r'<node id="(\d+)"', cblk)}
    sblk = re.search(r'<Surface name="curve_all">(.*?)</Surface>', t, re.S).group(1)
    inself = {tuple(sorted(int(v) for v in f.split(','))) for f in re.findall(r'>([^<]+)</quad4>', sblk)}
    ends = [f for f in ext if all(n in tied for n in f) and tuple(sorted(f)) not in inself]
    assert ends, 'no tied end faces found'
    surf = ('\t\t<Surface name="curve_endfaces">\n'
            + ''.join(f'\t\t\t<quad4 id="{k + 1}">{",".join(map(str, f))}</quad4>\n' for k, f in enumerate(ends))
            + '\t\t</Surface>\n')
    i = t.index('\t\t<Surface name="curve_all">')
    t = t[:i] + surf + t[i:]
    pair = ('\t\t<SurfacePair name="curve_skin_endfaces">\n\t\t\t<primary>curve_all</primary>\n'
            '\t\t\t<secondary>curve_endfaces</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name="curve_all_self">')
    t = t[:i] + pair + t[i:]
    c = re.search(r'\t\t<contact name="curve_all_self" surface_pair="curve_all_self" type="sliding-elastic">.*?</contact>\n', t, re.S).group(0)
    t = t.replace(c, c + c.replace('name="curve_all_self" surface_pair="curve_all_self"',
                                   'name="curve_skin_endfaces" surface_pair="curve_skin_endfaces"'))
    return t, [f'contact curve_skin_endfaces: curve skin (curve_all) against the {len(ends)} curve end faces tied into the '
               f'walls by canal_wrap_ends (new surface curve_endfaces), settings as curve_all_self (penalty 5, search '
               f'radius 0.5, one pass) [NOT IN SOURCE]']


if __name__ == '__main__':
    src, new = sys.argv[1:3]
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8', newline='').read()
    t, log = convert(t)
    print(new); print('\n'.join(log))
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batchC018.py)\nONE change (the mid-canal crossing):\n- '
          + '\n- '.join(log) + '\n')
