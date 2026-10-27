"""Batch C004 (2026-10-25, the colleague's turn): OUTPUT ONLY - must-point plot output, so a run's .xplt stays small.
L465-line runs write a plot state at every converged step (~2000-2700 steps x ~1 MB = 2-2.5 GB); FEBio 4.x ignored
plot_stride 5 (C002). Here: plot_level PLOT_MAJOR_ITRS -> PLOT_MUST_POINTS, and dtmax 0.005 -> a load curve that is 0.005
everywhere (so the maximum step is unchanged) with a point every DT (default 0.01) from 0 to the end time: the curve's points
become must points (skill febio-xml-format.md gotcha 31), so the solver lands exactly on t = 0.01, 0.02, ... and only those
states are written (~151 states, ~150 MB). No model change; the only solver effect is that a step is shortened when it would
step over a must point.
usage: py -3.10 build_batchC004.py SRC_RUN NEW_RUN [--dt 0.01]    (or import mp_output(t) in other builders)"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from paths import RUNS_DIR  # noqa: E402


def mp_output(t, dt=0.01):
    assert t.count('<dtmax>0.005</dtmax>') == 1 and t.count('<plot_level>PLOT_MAJOR_ITRS</plot_level>') == 1
    nsteps = int(re.search(r'<time_steps>(\d+)</time_steps>', t).group(1))
    step = float(re.search(r'<step_size>([^<]+)</step_size>', t).group(1))
    tend = nsteps * step
    lcid = max(int(v) for v in re.findall(r'<load_controller id="(\d+)"', t)) + 1
    t = t.replace('<dtmax>0.005</dtmax>', f'<dtmax lc="{lcid}">0.005</dtmax>')
    t = t.replace('<plot_level>PLOT_MAJOR_ITRS</plot_level>', '<plot_level>PLOT_MUST_POINTS</plot_level>')
    n = int(round(tend / dt))
    pts = ''.join(f'\t\t\t\t<pt>{k * dt:.4g},0.005</pt>\n' for k in range(n + 1))
    lc = (f'\t\t<load_controller id="{lcid}" name="dtmax_must_points" type="loadcurve">\n'
          '\t\t\t<interpolate>LINEAR</interpolate>\n\t\t\t<points>\n' + pts + '\t\t\t</points>\n\t\t</load_controller>\n')
    i = t.index('\t</LoadData>')
    t = t[:i] + lc + t[i:]
    return t, (f'output only: plot_level PLOT_MUST_POINTS, dtmax 0.005 -> load curve {lcid} "dtmax_must_points" (0.005 '
               f'everywhere, a point every {dt:g} to t {tend:g} = must points; ~{n + 1} plot states)')


if __name__ == '__main__':
    a = sys.argv[1:]
    dt = float(a[a.index('--dt') + 1]) if '--dt' in a else 0.01
    src, new = a[0], a[1]
    t = open(os.path.join(RUNS_DIR, src, src + '.feb'), encoding='utf-8', newline='').read()
    t, msg = mp_output(t, dt)
    print(new, msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batchC004.py)\nNO model change: {msg}\n')
