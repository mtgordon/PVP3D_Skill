"""Batch 26 (2026-09-24, user-chosen), on L11_stab (= the reference L9_fitall_edge without the 145 stab_* ground springs,
which are not in Abaqus; the user found its forward-moving vaginal wall realistic), fast base, to t = 1.0:
  L12_stabdamp     + the existing mass damping (settle_damping, C = 20/s) on from t = 0 instead of t = 1.0-1.05. NOT IN
                   SOURCE (numerical): without the springs the wall moves fast (up to 165-214 mm/s) and the runs crawl;
                   against a quasi-static ~15 mm/s motion C = 20/s is ~0.4 N, against > 130 N of applied pressure
  L12_stabdamp_p3  L12_stabdamp + every surface pressure at 1/3 of the source (user: "decrease the load on all surfaces to
                   match that of the levator"): Load-AVW, Load-PVW, Load-PeB-top, Load-top 0.014 -> 0.00467 MPa, like
                   Load-LA. NOT IN SOURCE
  L12_clusl_x1     + every CL/USL connector at full strength (x1: CL groups x0.2 and x0.1 and USL-bottom x0.1 given the
                   full-strength table) in the lofts' strip fits. NOT IN SOURCE: a sensitivity run, how much do the apex
                   and the LA depend on ligament strength now that the ground springs no longer hold the apex (they
                   carried ~3.0 N there at t = 1 against ~0.9 N for all four CL/USL)
usage: py -3.10 build_batch26.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

from variants6 import Model6, JOBS, emit
from loft_survey import survey
from loft_fit_all import fit_loft

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L11_stab', 'L11_stab.feb')
FULL = {'CL': 'ConnSect-ligment-CL-1-9-x%stiff', 'USL': 'ConnSect-ligment-USL-top-x%stiff'}
WANT = set(sys.argv[1:])


def want(name):
    return not WANT or name in WANT


def new(name):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(BASE)


def damping_from_start(m):
    lc = next(l for l in m.root.find('LoadData') if l.get('name') == 'settle_on')
    users = [e for e in m.root.iter() if e.get('lc') == lc.get('id')]
    assert [u.tag for u in users] == ['C'], [u.tag for u in users]
    pts = lc.find('points')
    old = [p.text for p in pts]
    for p in list(pts):
        pts.remove(p)
    for txt in ('0,1', '1.05,1'):
        ET.SubElement(pts, 'pt').text = txt
    m.log.append(f'NOT IN SOURCE (numerical): mass damping settle_damping (C = 20/s) on from t = 0: load curve settle_on '
                 f'{old} -> [0,1; 1.05,1] (it drives only the damping); Abaqus explicit has no such damping')


def pressures_third(m):
    for name in ('Load-AVW', 'Load-PVW', 'Load-PeB-top', 'Load-top'):
        m.scale_surface_load(name, 1 / 3)


def clusl_full_strength(m):
    f, parts, inst, asm, conns, behav, res, elem_of, bcs = survey(m.src)
    sprung = {frozenset(p) for pairs in f.discsets.values() for p in pairs}
    for fam in ('CL-L', 'CL-R', 'USL-L', 'USL-R'):
        loft, rows = res[fam]
        X = np.array([f.nodes[n] for n in sorted({n for c in f.elem_blocks[loft][1].values() for n in c})])
        on = [r for r in rows if np.linalg.norm(X - r['xb'], axis=1).min() < 0.01
              and frozenset((r['fa'], r['fb'])) not in sprung]
        full = FULL[fam.split('-')[0]]
        n_up = sum(r['behavior'] != full for r in on)
        r = fit_loft(f, loft, [dict(x, behavior=full) for x in on], behav, 0.49, 47.0)
        m.set_loft_ogden(loft, [(r['c'], r['m'])], 250.0, suffix='_x1', note=(
            f'NOT IN SOURCE (sensitivity): all {len(on)} connectors at full strength ({n_up} raised to {full}), strip fit '
            f"rms {100 * r['err']:.1f} % over u = 2-47 mm"))


if want('L12_stabdamp'):
    m = new('L12_stabdamp'); damping_from_start(m); emit('L12_stabdamp', m)
if want('L12_stabdamp_p3'):
    m = new('L12_stabdamp_p3'); damping_from_start(m); pressures_third(m); emit('L12_stabdamp_p3', m)
if want('L12_clusl_x1'):
    m = new('L12_clusl_x1'); clusl_full_strength(m); emit('L12_clusl_x1', m)
