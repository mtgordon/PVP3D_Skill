"""Batch 287 (2026-10-15; the user away a few hours: "keep trying things you think would be good").
The curves' fold-through at MAX is a SNAP: in L285_distal_off045_MAX one converged step (t 0.67709 -> 0.67868) moved a curve
surface node 2.0 mm and the penetrating nodes jumped 3 -> 42 (ordinary steps there move nodes 0.3-0.4 mm). No contact with a
0.5 mm search radius can catch a 2 mm jump. Here (a solver setting, not a model change): dtmax follows a load curve, 0.005
(as now) to t T0, DT from T0 + 0.01 to T1, 0.005 again from T1 + 0.01 (the curve's value IS dtmax and its points become must
points: tested on a 2-block file 10-15).
usage: py -3.10 build_batch287.py SRC_RUN NEW_RUN [--dt 0.00025] [--t0 0.60] [--t1 0.75]"""
import os
import re
import sys

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new, dt=0.00025, t0=0.60, t1=0.75):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    assert t.count('<dtmax>0.005</dtmax>') == 1
    lcid = max(int(v) for v in re.findall(r'<load_controller id="(\d+)"', t)) + 1
    t = t.replace('<dtmax>0.005</dtmax>', f'<dtmax lc="{lcid}">0.005</dtmax>')
    pts = [(0, 0.005), (t0, 0.005), (t0 + 0.01, dt), (t1, dt), (t1 + 0.01, 0.005), (1, 0.005)]
    lc = (f'\t\t<load_controller id="{lcid}" name="dtmax_fold_window" type="loadcurve">\n\t\t\t<interpolate>LINEAR</interpolate>\n'
          '\t\t\t<points>\n' + ''.join(f'\t\t\t\t<pt>{a},{b}</pt>\n' for a, b in pts) + '\t\t\t</points>\n\t\t</load_controller>\n')
    i = t.index('\t</LoadData>')
    t = t[:i] + lc + t[i:]
    d = os.path.join(RUNS, new)
    assert not os.path.exists(d), new
    os.makedirs(d)
    with open(os.path.join(d, new + '.feb'), 'w', encoding='utf-8', newline='') as f:
        f.write(t)
    with open(os.path.join(d, new + '.feb.changes.txt'), 'w', encoding='utf-8') as f:
        f.write(f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch287.py)\nsolver setting only: dtmax 0.005 -> load curve '
                f'{lcid} "dtmax_fold_window": {dt} from t {t0 + 0.01:.2f} to {t1:.2f} (0.005 outside; curve points are must points)\n')
    print(new, f'dtmax {dt} in t {t0 + 0.01:.2f}-{t1:.2f}, lc {lcid}')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, dflt: float(a[a.index(k) + 1]) if k in a else dflt   # noqa: E731
    main(a[0], a[1], g('--dt', 0.00025), g('--t0', 0.60), g('--t1', 0.75))
