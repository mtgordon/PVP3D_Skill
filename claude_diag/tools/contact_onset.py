"""Where a sliding contact makes a run crawl, for any contact pair (generalizes pvw_la_onset.py, 2026-09-27).

For each converged state in [T0, T1]: the solver's effort per step (log), how many facets of each side carry pressure,
and for the watched nodes (--nodes, else the fastest surface nodes) their speed, whether they reversed direction since
the previous step, and the pressure on the facets that hold them on each side (primary side = the pass that projects the
primary's points onto the secondary; secondary side = the two-pass second pass). Surfaces are found in the plot file by
name: <PAIR>Primary / <PAIR>Secondary (the converter's SlidingElastic1) or <PAIR>_primary / <PAIR>_secondary.
usage: py -3.10 contact_onset.py RUN PAIR T0 T1 [--nodes 2070,2071] [--every K] [--top N]
       PAIR may also be the two surfaces' names, primary first: SlidingElastic1Primary,canal_seam_strip
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt  # noqa: E402
from pvw_la_onset import log_steps  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('pair'); ap.add_argument('t0', type=float); ap.add_argument('t1', type=float)
    ap.add_argument('--nodes', default=''); ap.add_argument('--every', type=int, default=1)
    ap.add_argument('--top', type=int, default=6)
    a = ap.parse_args()
    x = Xplt(os.path.join(RUNS, a.run, a.run + '.xplt'))
    by = {s.get('name'): s for s in x.surfaces}
    sp = {}
    if ',' in a.pair:   # the two surfaces by name, primary first (a contact whose surfaces are not named after it)
        sp['primary'], sp['secondary'] = (by[n] for n in a.pair.split(','))
    for side, alts in (('primary', ('Primary', '_primary')), ('secondary', ('Secondary', '_secondary'))):
        if side not in sp:
            sp[side] = next(by[a.pair + s] for s in alts if a.pair + s in by)
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    nodes_of = {side: sorted({n for _, f in sp[side]['faces'] for n in f}) for side in sp}
    allsurf = sorted(set(nodes_of['primary']) | set(nodes_of['secondary']))
    watch = [idx[int(v)] for v in a.nodes.split(',') if v] if a.nodes else None
    facets_with = {side: {} for side in sp}
    for side in sp:
        for fi, (_, f) in enumerate(sp[side]['faces']):
            for n in f:
                facets_with[side].setdefault(n, []).append(fi)
    conv = [i for i, s in enumerate(x.states) if s[1] == 0] or list(range(len(x.states)))
    steps = log_steps(os.path.join(RUNS, a.run, a.run + '.log'))
    st_t = np.array([s[0] for s in steps])
    win = [k for k in conv if a.t0 <= x.states[k][0] <= a.t1][::a.every]
    print(f"{a.run} {a.pair}: primary {sp['primary']['nf']} facets, secondary {sp['secondary']['nf']} facets")
    prev_v = None
    for k in win:
        t = x.states[k][0]
        ci = conv.index(k)
        if ci == 0:
            continue
        kp = conv[ci - 1]
        u, up = x.var(k, 'displacement'), x.var(kp, 'displacement')
        v = (u - up) / (t - x.states[kp][0])
        spd = np.linalg.norm(v, axis=1)
        j = int(np.argmin(np.abs(st_t - t))) if len(st_t) else None
        dt, its, fails = (steps[j][1], steps[j][2], steps[j][3]) if j is not None and abs(st_t[j] - t) < 1e-5 else (0, '', '')
        cp = x.var(k, 'contact pressure') or {}
        pf = {}
        for side in sp:
            sid = sp[side]['id']
            pf[side] = (np.asarray(cp[sid]).reshape(sp[side]['nf'], -1).max(axis=1) if sid in cp
                        else np.zeros(sp[side]['nf']))
        w = watch if watch is not None else [allsurf[i] for i in np.argsort(-spd[allsurf])[:a.top]]
        cells = []
        for n in w:
            rev = ''
            if prev_v is not None:
                rev = 'R' if np.dot(v[n], prev_v[n]) < 0 else ''
            pp = max((pf['primary'][f] for f in facets_with['primary'].get(n, [])), default=np.nan)
            ps = max((pf['secondary'][f] for f in facets_with['secondary'].get(n, [])), default=np.nan)
            cells.append(f'{int(x.node_ids[n])}:{spd[n]:.0f}{rev} p{pp:.3g}/s{ps:.3g}')
        print(f"t {t:.5f} dt {dt:.1e} its {its} fails {fails} | touching p {int((pf['primary'] > 1e-9).sum())} "
              f"s {int((pf['secondary'] > 1e-9).sum())} | vmax surf {spd[allsurf].max():.0f} mm/s | " + '  '.join(cells))
        prev_v = v


if __name__ == '__main__':
    main()
