"""Batch C003 (2026-10-24, the colleague's turn; Mark's idea 1 for the corner, NEXT_SESSION_PROMPT_2026-10-24): a local
stiffer band in the curve elements that fold. The curve self-crossings sit at rest z -31 (left) / -21 (right); every curve
hex (canal_wrap_hex_avw / _pvw / _pvw_fine; not the cervix part) whose rest centroid has ZLO <= z < ZHI moves to a twin
domain <parent>_band with Curve_band = Curve_ogmid with c1, c2 and k x SCALE (same Ogden exponents, so the same
tension / compression ratio, just stiffer: more bending stiffness, a wider fold radius). The rest of the curves and the
whole wall keep the soft middle law. Element ids, nodes, surfaces and contacts are unchanged; each band domain gets its
parent's mass damping. The domain names still start with canal_wrap_hex, so curve_cross.py & co. see them.
NOT IN SOURCE (the source has no curves).
usage: py -3.10 build_batchC003.py SRC_RUN NEW_RUN [--zlo -33] [--zhi -19] [--scale 2]"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from paths import RUNS_DIR  # noqa: E402

PARENTS = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')


def convert(t, zlo, zhi, scale, tag='band'):
    X = {int(i): np.array([float(v) for v in xyz.split(',')])
         for i, xyz in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    log, counts = [], {}
    for p in PARENTS:
        m = re.search(r'(\t\t<Elements type="hex8" name="%s">\n)(.*?)(\t\t</Elements>\n)' % p, t, re.S)
        keep, band = [], []
        for line in m.group(2).splitlines(keepends=True):
            ids = [int(v) for v in re.search(r'>([^<]+)<', line).group(1).split(',')]
            z = np.mean([X[i][2] for i in ids])
            (band if zlo <= z < zhi else keep).append(line)
        counts[p] = len(band)
        if not band:
            continue
        if not keep:   # the whole domain is in the band: switch its material (an empty <Elements> block fails in FEBio)
            t = t.replace('\t\t<SolidDomain name="%s" mat="Curve_ogmid" />\n' % p,
                          '\t\t<SolidDomain name="%s" mat="Curve_%s" />\n' % (p, tag))
            continue
        new = (m.group(1) + ''.join(keep) + m.group(3)
               + '\t\t<Elements type="hex8" name="%s_%s">\n' % (p, tag) + ''.join(band) + m.group(3))
        t = t.replace(m.group(0), new)
        t = t.replace('\t\t<SolidDomain name="%s" mat="Curve_ogmid" />\n' % p,
                      '\t\t<SolidDomain name="%s" mat="Curve_ogmid" />\n'
                      '\t\t<SolidDomain name="%s_%s" mat="Curve_%s" />\n' % (p, p, tag, tag))
        b = re.search(r'\t\t<body_load name="wall_damping_%s" type="mass damping" elem_set="%s">.*?</body_load>\n' % (p, p), t, re.S)
        t = t.replace(b.group(0), b.group(0) + b.group(0).replace(p, p + '_' + tag))
    assert 1 <= t.count('mat="Curve_%s"' % tag) <= 3
    m = re.search(r'(\t\t<material id="(\d+)" name="Curve_ogmid" type="Ogden">.*?</material>\n)', t, re.S)
    blk = m.group(1)
    ids = [int(v) for v in re.findall(r'<material id="(\d+)"', t)]
    nb = blk.replace('id="%s"' % m.group(2), 'id="%d"' % (max(ids) + 1)).replace('Curve_ogmid', 'Curve_' + tag)
    for k in ('c1', 'c2', 'k'):
        old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), nb).group(1))
        nb = re.sub(r'<%s>[^<]+</%s>' % (k, k), '<%s>%.7g</%s>' % (k, old * scale, k), nb)
    i = t.index('\t</Material>')
    t = t[:i] + nb + t[i:]
    log.append('curve band: %s curve hexes with rest centroid %g <= z < %g (%s) -> twin domains *_%s with Curve_%s = '
               'Curve_ogmid c1 / c2 / k x%g (same m1 / m2), parent damping copied [NOT IN SOURCE]'
               % (sum(counts.values()), zlo, zhi, ', '.join('%s %d' % kv for kv in counts.items()), tag, tag, scale))
    return t, log


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('new')
    ap.add_argument('--zlo', type=float, default=-33); ap.add_argument('--zhi', type=float, default=-19)
    ap.add_argument('--scale', type=float, default=2); ap.add_argument('--tag', default='band')
    a = ap.parse_args()
    t = open(os.path.join(RUNS_DIR, a.src, a.src + '.feb'), encoding='utf-8', newline='').read()
    t, log = convert(t, a.zlo, a.zhi, a.scale, a.tag)
    print(a.new); print('\n'.join(log))
    write(a.new, t, f'base: claude_diag/runs/{a.src}/{a.src}.feb (tools/build_batchC003.py)\nONE change (Mark\'s corner idea 1, a '
                    'local stiffer band where the curve folds):\n- ' + '\n- '.join(log) + '\n')
