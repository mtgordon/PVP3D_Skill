"""Comparison table for batch 18 (every loft against its own connectors): convergence, LA displacement, and each
connector family's end elongation (the in-situ check of a loft against its Abaqus connectors).

usage: py -3.10 batch18_table.py RUN [RUN ...]
Per run: last converged t, converged steps, failed attempts, negative-jacobian warnings, status, wall time (log
'Total elapsed time', else log creation -> last write); then at t = 1.0 and at the last converged state: LA |u| median /
p90 / max, max |u| of the rest, the fastest nodal speed; and per connector family the median / max elongation between
the live end nodes, with the source force it implies (sum over the family, each connector's own table).
"""
import os
import re
import sys
import time

import numpy as np

import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR
from xplt import Xplt
from febmodel import Feb
from loft_survey import survey

FAMS = ('AVW-Para-L', 'AVW-Para-R', 'CL-L', 'CL-R', 'USL-L', 'USL-R', 'PM', 'PM_PeB_Left_', 'PeB-constrin',
        'P-arcus-L', 'P-arcus-R')


def log_info(run):
    p = os.path.join(RUNS_DIR, run, run + '.log')
    s = open(p, encoding='latin-1').read()
    conv = re.findall(r'^------- converged at time : ([0-9.eE+-]+)', s, re.M)
    fails = len(re.findall(r'failed to converge at time', s))
    negj = len(re.findall(r'negative jacobian', s))
    st = 'NORMAL' if 'N O R M A L   T E R M' in s else ('ERROR' if 'E R R O R   T E R M' in s else 'running')
    m = re.search(r'Total elapsed time \.+ : (\S+)', s)
    if m:
        wall = m.group(1)
    else:
        stt = os.stat(p)
        sec = stt.st_mtime - stt.st_ctime
        wall = time.strftime('%H:%M:%S', time.gmtime(sec)) + ('' if st == 'running' else ' (no summary)')
    return float(conv[-1]) if conv else 0.0, len(conv), fails, negj, st, wall


def main(runs):
    base_feb = os.path.join(RUNS_DIR, runs[0], runs[0] + '.feb')
    f0, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(base_feb)
    tabs = {b: np.array(behav[b][0]['table']) for b in behav if behav[b]}
    print('| run | t reached | steps | failed attempts | neg. J | status | wall |')
    print('|---|---|---|---|---|---|---|')
    info = {}
    for r in runs:
        info[r] = log_info(r)
        t, n, fl, nj, st, wall = info[r]
        print(f'| {r} | {t:.4f} | {n} | {fl} | {nj} | {st} | {wall} |')
    print()
    for r in runs:
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
        la = set()
        for i, d in enumerate(x.domains):
            if i < len(names) and names[i].startswith('LA_'):
                la.update(int(nd) for c in d['conn'] for nd in c)
        la = np.array(sorted(la))
        rest = np.setdiff1d(np.arange(len(x.node_ids)), la)
        conv = [i for i, s in enumerate(x.states) if s[1] == 0]
        times = np.array([x.states[i][0] for i in conv])
        picks = sorted({int(np.argmin(abs(times - 1.0)))} | {len(times) - 1}) if times[-1] >= 0.99 else [len(times) - 1]
        for k in picks:
            u = x.var(conv[k], 'displacement')
            un = np.linalg.norm(u, axis=1)
            k0 = max(k - 1, 0)
            v = np.linalg.norm(u - x.var(conv[k0], 'displacement'), axis=1) / max(times[k] - times[k0], 1e-12)
            q = np.percentile(un[la], [50, 90])
            line = (f'{r:12s} t {times[k]:.3f}: LA {q[0]:5.1f} / {q[1]:5.1f} / {un[la].max():5.1f} mm, rest max '
                    f'{un[rest].max():5.1f} mm, max speed {v.max():7.1f} mm/s | elongation med/max [mm] (source force N):')
            parts_ = []
            for fam in FAMS:
                if fam not in res:
                    continue
                e = []
                F = 0.0
                for row in res[fam][1]:
                    pts = []
                    for xyz in (row['xa'], row['xb']):
                        cand = ids[np.linalg.norm(P - xyz, axis=1) < 1e-6]
                        lv = [nd for nd in cand if nd in live]
                        nd = lv[0] if lv else cand[0]
                        pts.append(xyz + u[idx[nd]])
                    el = np.linalg.norm(pts[1] - pts[0]) - row['L']
                    e.append(el)
                    tb = tabs[row['behavior']]
                    F += np.interp(el, tb[:, 1], tb[:, 0]) if el > 0 else 0.0
                e = np.array(e)
                parts_.append(f'{fam} {np.median(e):.1f}/{e.max():.1f} ({F:.3g})')
            print(line + ' ' + '; '.join(parts_))


if __name__ == '__main__':
    main(sys.argv[1:])
