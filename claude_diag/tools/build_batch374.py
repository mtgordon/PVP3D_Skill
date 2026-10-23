"""Batch 374 (2026-10-21; the user: "keeping the original material strength up to about 35 and then making it get very stiff very
quickly ... do the same for CL. Let's go even steeper than the healthy slope for now.").
NOT IN SOURCE (a new connector law; LOCAL connectors, CL + USL): every CL / USL spring table becomes, in absolute force (scale 1):
  - up to X0 (35 mm): the source table x SCALE (the case's impairment, MAX-lite 0.2), interpolated at X0;
  - X0 .. X0 + RAMP (4 mm): the slope rises linearly from the impaired source slope to S1 (a smooth quadratic, no kink);
  - beyond: straight at S1 to a far point (200 mm) with extend constant (no plateau in range).
S1 = SF x the HEALTHY table's end slope of that bundle ((F(51.9) - F(48.4)) / 3.5 at x1: 0.254 N/mm for CL 1-9 / USL top,
/5 and /10 for the weaker tables), so the three source tables keep their ratios. Tension only (F = 0 at <= 0) as before.
usage: py -3.10 build_batch374.py SRC_RUN NEW_RUN [--x0 35] [--ramp 4] [--sf 3] [--scale 0.2]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
MATS = [f'{a}-{s}_conn_{i}_mat' for a, n in (('CL', 3), ('USL', 2)) for s in 'LR' for i in range(1, n + 1)]


def law(xs, fs, scale, x0, ramp, sf, far=200.0):
    xs, fs = np.array(xs), np.array(fs)
    f0 = scale * np.interp(x0, xs, fs)
    j = np.searchsorted(xs, x0)                       # segment holding x0
    s0 = scale * (fs[j] - fs[j - 1]) / (xs[j] - xs[j - 1])
    s1 = sf * (fs[-1] - fs[-2]) / (xs[-1] - xs[-2])   # x sf the healthy end slope
    pts = [(x, scale * f) for x, f in zip(xs, fs) if x < x0] + [(x0, f0)]
    for d in np.arange(1, ramp + 1):
        pts.append((x0 + d, f0 + s0 * d + (s1 - s0) * d * d / (2 * ramp)))
    xe, fe = pts[-1]
    pts.append((far, fe + s1 * (far - xe)))
    return pts, s0, s1


def main(src, new, x0, ramp, sf, scale):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    log = []
    for m in MATS:
        mm = re.search(rf'(<discrete_material id="\d+" name="{re.escape(m)}" type="nonlinear spring">)(.*?)(</discrete_material>)', t, re.S)
        body = mm.group(2)
        assert abs(float(re.search(r'<scale>([^<]+)</scale>', body).group(1)) - scale) < 1e-9, m
        P = [tuple(float(v) for v in p.split(',')) for p in re.findall(r'<pt>([^<]+)</pt>', body)]
        pts, s0, s1 = law([p[0] for p in P], [p[1] for p in P], scale, x0, ramp, sf)
        rows = ''.join(f'\t\t\t\t\t<pt>{x:g},{f:.6g}</pt>\n' for x, f in pts)
        nb = re.sub(r'<scale>[^<]+</scale>', '<scale>1</scale>', body)
        nb = re.sub(r'(<points>\n).*?(\t*</points>)', lambda q: q.group(1) + rows + q.group(2), nb, flags=re.S)
        t = t[:mm.start(2)] + nb + t[mm.end(2):]
        log.append(f'{m}: F({x0:g}) {pts[len([p for p in P if p[0] < x0])][1]:.4g} N, slope {s0:.4g} -> {s1:.4g} N/mm')
    msg = (f'NOT IN SOURCE (new CL + USL law, LOCAL connectors): source table x{scale:g} up to {x0:g} mm, then the slope ramps over '
           f'{ramp:g} mm to {sf:g} x the healthy table\'s end slope and stays straight (no plateau); scale 1 with absolute forces '
           f'(tools/build_batch374.py): ' + '; '.join(log))
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch374.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k_, d: float(a[a.index(k_) + 1]) if k_ in a else d   # noqa: E731
    main(a[0], a[1], g('--x0', 35.0), g('--ramp', 4.0), g('--sf', 3.0), g('--scale', 0.2))
