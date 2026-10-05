"""Batch 31 (2026-09-25, overnight): the all-connector model, which passes the source-pressure wall.

Batches 29-30 (the stab-spring-free model at the SOURCE pressures, fast base, 1 s): every model with lofts crawled at
t ~ 0.65-0.84 (39-523 failed attempts); L16_allconn2 (every Abaqus CONN3D2 connector as a FEBio spring, no lofts: the
source's supports exactly) was at t 0.80 with 1 failed attempt after 25 min. One change each:
  L17D_allconn      L16_allconn2 like for like: the arcus chain mass at the Abaqus truss density (chain_mass_mat 0.0011 ->
                    0.00011: the source's 33 kg instead of the fast base's 329 kg). With the source's 1 s ramp, source
                    pressures, no stab springs, no damping, velocity plotted: the candidate for the Abaqus comparison
  L17_lofts4        L16_nzconn (the five near-zero lofts as connectors) + the P-arcus lofts as their connectors
                    (+ the chain unsplit): keeps the user's four fitted lofts AVW-Para L/R, CL L/R, USL L/R and PM. Are
                    the soft lofts (near-zero, P-arcus) the whole problem?
  L17_allconn_hold  L16_allconn2 run on to t = 2.5 with the load held; the settle damping (on at t = 1.0-1.05) at
                    C = 2/s instead of 20 (NOT IN SOURCE, after the ramp only): the equilibrium at the source loads.
                    L14_p3_hold showed C = 20/s on the fast base's 329 kg chain makes the approach a slow creep
                    (C x chain mass ~ 30 N per 5 mm/s)
usage: py -3.10 build_batch31.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

NAMES = ('L17D_allconn', 'L17_lofts4', 'L17_allconn_hold')
if not sys.argv[1:]:
    sys.argv += list(NAMES)   # build_batch30 reads the same sys.argv on import and then builds none of its own runs

from variants6 import Model6, JOBS, emit  # noqa: E402
from build_batch30 import parcus_to_connectors, unsplit_chain  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
ALLCONN = os.path.join(RUNS, 'L16_allconn2', 'L16_allconn2.feb')
NZCONN = os.path.join(RUNS, 'L16_nzconn', 'L16_nzconn.feb')
WANT = set(sys.argv[1:])


def want(name):
    return not WANT or name in WANT


def fresh(name, base):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(base)


def chain_mass_source(m):
    mat = next(x for x in m.root.find('Material') if x.get('name') == 'chain_mass_mat')
    old = mat.find('density').text
    mat.find('density').text = '0.00011'
    m.log.append(f'chain_mass_mat density {old} -> 0.00011: the Abaqus truss density (ATLA-Hyper, Posterior_Arcus-Hyper, '
                 f'area 1 mm2), like for like; the fast base had x10')


def settle_c(m, c=2.0):
    bl = next(b for b in m.root.find('Loads') if b.get('name') == 'settle_damping')
    el = bl.find('C')
    old = el.text
    el.text = '%g' % c
    lc = next(l for l in m.root.find('LoadData') if l.get('id') == el.get('lc'))
    m.log.append(f'NOT IN SOURCE (after the ramp only): settle_damping C {old} -> {c:g}/s, switched on by '
                 f'{lc.get("name")} {[p.text for p in lc.find("points")]}')


if want('L17D_allconn'):
    m = fresh('L17D_allconn', ALLCONN); chain_mass_source(m); emit('L17D_allconn', m)
if want('L17_lofts4'):
    m = fresh('L17_lofts4', NZCONN); parcus_to_connectors(m); unsplit_chain(m); emit('L17_lofts4', m)
if want('L17_allconn_hold'):
    m = fresh('L17_allconn_hold', ALLCONN); m.set_end_time(2.5)
    m.log.append('the load is held from t = 1 to 2.5 (Amp-1 extends CONSTANT)'); settle_c(m, 2.0)
    emit('L17_allconn_hold', m)


# Batch 32 (2026-09-25 ~01:00): the next wall is the LA_ICM snap-through near t ~ 0.9 (97 % load; L16_allconn2: LA
# 28.3 / 38.5 / 42.1 mm, LA_ICM nodes at 1012 mm/s). The real cutback (aggressiveness 1, cutback 0.5) cut the failed
# attempts 2-5x in every crawl so far (L16_aggr vs L15_src, L16_conn_aggr vs L15_conn) and did nothing in smooth runs
# (L8_aggr), so it goes on the connector models: solver only, one change each.
#   L18_allconn_aggr   L16_allconn2 + the cutback
#   L18_lofts4_aggr    L17_lofts4 + the cutback
#   L18D_allconn_aggr  L17D_allconn (like for like) + the cutback
#   L18_hold_aggr      L17_allconn_hold (hold to t = 2.5, settle C = 2/s after t = 1) + the cutback
def cutback(m):
    m.set_solver(aggressiveness=1, cutback=0.5)
    m.log.append('solver only: a real cutback (FEBio reads cutback only with aggressiveness 1): each retry halves the '
                 'step instead of taking dt0/21 off it')


for _name, _base in (('L18_allconn_aggr', 'L16_allconn2'), ('L18_lofts4_aggr', 'L17_lofts4'),
                     ('L18D_allconn_aggr', 'L17D_allconn'), ('L18_hold_aggr', 'L17_allconn_hold')):
    if _name in WANT:
        m = fresh(_name, os.path.join(RUNS, _base, _base + '.feb')); cutback(m); emit(_name, m)

# L18D_lofts4_aggr (2026-09-25 ~01:30): L18_lofts4_aggr like for like (chain mass at the Abaqus density): the user's four
# fitted lofts on the Abaqus-comparison base. One change from L18_lofts4_aggr.
if 'L18D_lofts4_aggr' in WANT:
    m = fresh('L18D_lofts4_aggr', os.path.join(RUNS, 'L18_lofts4_aggr', 'L18_lofts4_aggr.feb')); chain_mass_source(m)
    emit('L18D_lofts4_aggr', m)
