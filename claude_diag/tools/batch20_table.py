"""Comparison table for batch 20 (solver settings from the 3-surface test, near-zero loft fits, tie penalty).

usage: py -3.10 batch20_table.py [--t T] RUN [RUN ...]   (--t: also report at the converged state nearest T, to
                                                          compare runs that stopped at different times)
1. Per run: change, last converged t, converged steps, failed attempts, negative-jacobian warnings, status, wall time.
2. At t = 1.0 (when reached) and at the last converged state: LA |u| median / p90 / max, and per connector family
   the median / max elongation between the live end nodes with the source force it implies (as batch18_table.py,
   plus PM_PeB right and PM_avw_bottom).
3. Where a run that did not finish was failing: at its last converged state, the domains of the 50 fastest nodes
   (speed over the last converged step) and the three domains with the lowest relative volume J.
4. The near-zero lofts' J (min / median / max): k = 250 mu0 was chosen for ~3x stretch, their short strips reach ~4-12x.
"""
import os
import sys
from collections import Counter, defaultdict

import numpy as np

import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR
from xplt import Xplt
from febmodel import Feb
from loft_survey import survey
from batch18_table import log_info

CHANGES = {
    'L6_ctl': 'control: L5X_all to t = 1',
    'L6_ls': 'line search on (lstol 0 -> 0.9)',
    'L6_fn': 'full Newton (max_ups 10 -> 0)',
    'L6_blk': "3-surface test's solver block",
    'L6_pmpeb': 'PM_PeB + PM_avw_bottom lofts fitted (near zero), tissue density',
    'L6_pebc': 'PeB-constrin loft fitted (near zero), tissue density',
    'L6_fitall': 'every non-P-arcus loft fitted',
    'L6_tie1': 'tnof tie penalty 100 -> 1 N/mm',
    'L6D_x': 'like-for-like (DM1) + LA load x 1/3 + the three fixes',
}
FAMS = ('AVW-Para-L', 'AVW-Para-R', 'CL-L', 'CL-R', 'USL-L', 'USL-R', 'PM', 'PM_PeB_Left_', 'PM_PeB_Right_',
        'PM_avw_bottom_left_', 'PM_avw_bottom_right_', 'PeB-constrin', 'P-arcus-L', 'P-arcus-R')
NEAR_ZERO = ('PM_PeB_Left_fan', 'PM_PeB_Right_fan', 'PM_avw_bottom_left_fan', 'PM_avw_bottom_right_fan',
             'PeB-constrin_fan')


def main(runs, t_at=None):
    f0, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(os.path.join(RUNS_DIR, runs[0], runs[0] + '.feb'))
    tabs = {b: np.array(behav[b][0]['table']) for b in behav if behav[b]}
    print('| run | change | t reached | steps | failed attempts | neg. J | status | wall |')
    print('|---|---|---|---|---|---|---|---|')
    info = {r: log_info(r) for r in runs}
    for r in runs:
        t, n, fl, nj, st, wall = info[r]
        print(f'| {r} | {CHANGES.get(r, "")} | {t:.4f} | {n} | {fl} | {nj} | {st} | {wall} |')
    for r in runs:
        print(f'\n== {r}')
        fr = Feb(os.path.join(RUNS_DIR, r, r + '.feb'))
        names = [d.get('name') for d in fr.root.find('MeshDomains')]
        live = set()
        for nm, (et, d) in fr.elem_blocks.items():
            for c in d.values():
                live.update(c)
        for pairs in fr.discsets.values():
            for a, b in pairs:
                live.update((a, b))
        ids = np.array(sorted(fr.nodes))
        P = np.array([fr.nodes[i] for i in ids])
        x = Xplt(os.path.join(RUNS_DIR, r, r + '.xplt'))
        idx = {int(nd): i for i, nd in enumerate(x.node_ids)}
        node_dom = defaultdict(set)      # plot node index (the plot's connectivity is 0-based indices) -> domains
        la = set()
        for i, d in enumerate(x.domains):
            nm = names[i] if i < len(names) else f'region{i + 1}'
            for j in np.unique(d['conn'].ravel()):
                node_dom[int(j)].add(nm)
            if nm.startswith('LA_'):
                la.update(int(j) for c in d['conn'] for j in c)
        la = np.array(sorted(la))
        conv = [i for i, s in enumerate(x.states) if s[1] == 0]
        if not conv:
            print('   no converged state in the plot file')
            continue
        times = np.array([x.states[i][0] for i in conv])
        picks = {int(np.argmin(abs(times - 1.0)))} | {len(times) - 1} if times[-1] >= 0.99 else {len(times) - 1}
        if t_at is not None:
            picks.add(int(np.argmin(abs(times - t_at))))
        picks = sorted(picks)
        for k in picks:
            u = x.var(conv[k], 'displacement')
            un = np.linalg.norm(u, axis=1)
            q = np.percentile(un[la], [50, 90])
            print(f'   t {times[k]:.3f}: LA {q[0]:.1f} / {q[1]:.1f} / {un[la].max():.1f} mm (median / p90 / max)')
            parts_ = []
            for fam in FAMS:
                if fam not in res:
                    continue
                e, F = [], 0.0
                for row in res[fam][1]:
                    pts = []
                    for xyz in (row['xa'], row['xb']):
                        cand = ids[np.linalg.norm(P - xyz, axis=1) < 1e-6]
                        lv = [nd for nd in cand if nd in live]
                        nd = lv[0] if lv else cand[0]
                        pts.append(xyz + (u[idx[nd]] if nd in idx else 0.0))
                    el = np.linalg.norm(pts[1] - pts[0]) - row['L']
                    e.append(el)
                    tb = tabs[row['behavior']]
                    F += np.interp(el, tb[:, 1], tb[:, 0]) if el > 0 else 0.0
                e = np.array(e)
                parts_.append(f'{fam} {np.median(e):.1f}/{e.max():.1f} ({F:.3g} N)')
            print('   connector-end elongation median/max mm (source force): ' + '; '.join(parts_))
            J = x.var(conv[k], 'relative volume')
            nz = []
            for i, nm in enumerate(names):
                if nm in NEAR_ZERO and J is not None and (i + 1) in J:
                    a = J[i + 1]
                    nz.append(f'{nm} {a.min():.3f}/{np.median(a):.3f}/{a.max():.3f}')
            print('   near-zero lofts J min/median/max: ' + '; '.join(nz))
        if info[r][4] != 'NORMAL' and len(conv) > 1:
            i1, i0 = conv[-1], conv[-2]
            u1, u0 = x.var(i1, 'displacement'), x.var(i0, 'displacement')
            v = np.linalg.norm(u1 - u0, axis=1) / max(times[-1] - times[-2], 1e-12)
            top = np.argsort(v)[::-1][:50]
            doms = Counter(nm for j in top for nm in node_dom[int(j)])
            J = x.var(i1, 'relative volume')
            low = sorted((float(J[i + 1].min()), nm) for i, nm in enumerate(names) if J is not None and (i + 1) in J)[:3]
            print(f'   failing at t {times[-1]:.4f}: max speed {v.max():.1f} mm/s; 50 fastest nodes in '
                  + ', '.join(f'{k} {c}' for k, c in doms.most_common(4))
                  + '; lowest J: ' + ', '.join(f'{nm} {j:.3f}' for j, nm in low))


if __name__ == '__main__':
    a = sys.argv[1:]
    t_at = None
    if a and a[0] == '--t':
        t_at, a = float(a[1]), a[2:]
    main(a, t_at)
