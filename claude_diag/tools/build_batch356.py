"""Batch 356 (2026-10-20; the user, on the curves collapsing on themselves at t ~0.83 in L352: option 3 = a higher bulk
modulus on the curves' bend elements only). L286's pictures (10-15) showed the elements on the inside of the fold losing volume
(V/V0 to 0.41): the curves' Yeoh k 0.21 is small against the c2-stiffened shear at large strain. L273 (tube line, k x5 on the
distal 15 mm) still folded, so a larger factor here.
NOT IN SOURCE (a local material change, LOCAL band): every hex of DOMS whose centre is at rest z < ZCUT moves to a new domain
<dom>_distal with material Curve_soft_distalk = Curve_soft with k x KF (c1, c2, density unchanged); the new domains get the same
mass damping as their parent domain. Surfaces, ties and contacts are node-based and stay as they are.
usage: py -3.10 build_batch356.py SRC_RUN NEW_RUN [--zcut -15] [--kf 20] [--df 1]
--df (2026-10-21): density x DF as well (the jitter recipe: local mass under the parent's mass damping; L366 = --kf 1 --df 100)"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
DOMS = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')


def main(src, new, zcut=-15.0, kf=20.0, df=1.0):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    m = re.search(r'\t\t<material id="(\d+)" name="Curve_soft" type="Yeoh">.*?</material>\n', t, re.S)
    k0 = float(re.search(r'<k>([^<]+)</k>', m.group(0)).group(1))
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    rho0 = float(re.search(r'<density>([^<]+)</density>', m.group(0)).group(1))
    newmat = re.sub(r'<k>[^<]+</k>', f'<k>{k0 * kf:g}</k>', m.group(0))
    newmat = re.sub(r'<density>[^<]+</density>', f'<density>{rho0 * df:g}</density>', newmat).replace(f'id="{m.group(1)}" name="Curve_soft"', f'id="{mid}" name="Curve_soft_distalk"')
    i = t.index('\t</Material>'); t = t[:i] + newmat + t[i:]   # at the end, in id order
    log = []
    for d in DOMS:
        blk = re.search(rf'(\t\t<Elements type="hex8" name="{d}">\n)(.*?)(\t\t</Elements>\n)', t, re.S)
        rows = re.findall(r'\t\t\t<elem id="\d+">[0-9,]+</elem>\n', blk.group(2))
        dist = [r for r in rows if np.mean([X[int(v)][2] for v in re.search(r'>([0-9,]+)<', r).group(1).split(',')]) < zcut]
        if not dist:
            log.append(f'{d}: 0'); continue
        keep = [r for r in rows if r not in set(dist)]
        if not keep:   # the whole domain is in the zone (canal_wrap_hex_pvw_fine): switch its material, no split
            dom = re.search(rf'\t\t<SolidDomain name="{d}" mat="Curve_soft" />\n', t).group(0)
            t = t.replace(dom, dom.replace('mat="Curve_soft"', 'mat="Curve_soft_distalk"'), 1)
            log.append(f'{d}: all {len(rows)} (domain switched)'); continue
        t = t.replace(blk.group(0), blk.group(1) + ''.join(keep) + blk.group(3) +
                      blk.group(1).replace(f'name="{d}"', f'name="{d}_distal"') + ''.join(dist) + blk.group(3), 1)
        dom = re.search(rf'\t\t<SolidDomain name="{d}" mat="Curve_soft" />\n', t).group(0)
        t = t.replace(dom, dom + f'\t\t<SolidDomain name="{d}_distal" mat="Curve_soft_distalk" />\n', 1)
        bl = re.search(rf'\t\t<body_load name="wall_damping_{d}" type="mass damping" elem_set="{d}">.*?</body_load>\n', t, re.S).group(0)
        t = t.replace(bl, bl + bl.replace(f'wall_damping_{d}"', f'wall_damping_{d}_distal"').replace(f'elem_set="{d}"', f'elem_set="{d}_distal"'), 1)
        log.append(f'{d}: {len(dist)} of {len(rows)}')
    msg = (f'NOT IN SOURCE (local band): the curves\' bend elements (hex centre at rest z < {zcut:g}) on Curve_soft_distalk = Curve_soft '
           f'with k x{kf:g} ({k0:g} -> {k0 * kf:g}), density x{df:g} ({rho0:g} -> {rho0 * df:g}) (c1, c2 unchanged) in new domains *_distal with the parent\'s mass damping '
           f'(tools/build_batch356.py): ' + '; '.join(log))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch356.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: float(a[a.index(k_) + 1]) if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--zcut', -15.0), g('--kf', 20.0), g('--df', 1.0))
