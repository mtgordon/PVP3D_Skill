"""Batch 33 (2026-09-25 ~02:25, overnight): the PVW-LA contact onset, the wall after the chain snap on the fast base.

L18_allconn_aggr (every Abaqus connector a spring, no lofts, real cutback, source pressures, no stab springs) passed the
chain snap (t ~ 0.9) and then crawled at t 0.907 with all 50 fastest nodes in the posterior vaginal wall (_PickedSet64)
at up to 6 m/s, steps down to 6e-5: from t 0.9035 ONE facet of the PVW_LA sliding-elastic contact (penalty 5,
auto_penalty, two_pass, offset 2) is in contact (p 0.012-0.016 MPa), none before. The four-loft model crawls where its
PVW_LA contact starts too (1-2 facets from t 0.95, p up to 0.09 MPa). One change each on L18_allconn_aggr:
  L19_pen5       PVW_LA penalty 5 -> 0.5 (auto_penalty stays on: a 10x softer contact onset)
  L19_nopvwla    PVW_LA contact removed. NOT IN SOURCE (Abaqus has this contact pair): a diagnostic only, does the run go
                 through without the onset?
usage: py -3.10 build_batch33.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L18_allconn_aggr', 'L18_allconn_aggr.feb')
WANT = set(sys.argv[1:])


def want(name):
    return not WANT or name in WANT


def fresh(name):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(BASE)


if want('L19_pen5'):
    m = fresh('L19_pen5')
    assert m.set_contact('PVW_LA', penalty=0.5) == 1
    m.log.append('PVW_LA: penalty 5 -> 0.5 (auto_penalty 1 kept: 10x softer); the source contact pair is kept')
    emit('L19_pen5', m)
if want('L19_nopvwla'):
    m = fresh('L19_nopvwla')
    con = m.root.find('Contact')
    c = next(c for c in con if c.get('name') == 'PVW_LA')
    con.remove(c)
    mesh = m.root.find('Mesh')
    sp = next(s for s in mesh.findall('SurfacePair') if s.get('name') == c.get('surface_pair'))
    mesh.remove(sp)
    m.log.append('NOT IN SOURCE (diagnostic): the PVW_LA sliding-elastic contact and its surface pair removed (the '
                 'surfaces stay, unused); Abaqus has this contact pair')
    emit('L19_nopvwla', m)

# L20_lofts4_hold (2026-09-25 ~02:50): the equilibrium at the source loads of the first model to reach them
# (L18_lofts4_aggr: t = 1, normal termination, LA 33.3 / 46.7 / 51.9 mm): held to t = 2.5, the settle damping (on at
# t = 1.0-1.05) at C = 2/s instead of 20 (NOT IN SOURCE, after the ramp only). L14_p3_hold: C = 20/s on the fast base's
# 329 kg chain made the approach a creep (tau 6-20 s, so a mode near omega ~ 1.4 rad/s; C ~ 2 omega = 2.8/s is critical).
if 'L20_lofts4_hold' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L20_lofts4_hold')), 'L20_lofts4_hold exists; not overwriting'
    m = Model6(os.path.join(RUNS, 'L18_lofts4_aggr', 'L18_lofts4_aggr.feb'))
    m.set_end_time(2.5)
    m.log.append('the load is held from t = 1 to 2.5 (Amp-1 extends CONSTANT)')
    bl = next(b for b in m.root.find('Loads') if b.get('name') == 'settle_damping')
    lc = next(l for l in m.root.find('LoadData') if l.get('id') == bl.find('C').get('lc'))
    assert [p.text for p in lc.find('points')] == ['1.0,0', '1.05,1'], [p.text for p in lc.find('points')]
    old = bl.find('C').text
    bl.find('C').text = '2'
    m.log.append(f'NOT IN SOURCE (after the ramp only): settle_damping C {old} -> 2/s, on at t = 1.0-1.05 (settle_on)')
    emit('L20_lofts4_hold', m)

# L21_lofts4_pen (2026-09-25 ~03:15): L19_nopvwla (no PVW_LA contact) reached t = 1 in 46 min with 3 failed attempts where
# its control crawled (the contact touches at 1-6 nodes near the end: without it the wall comes within 0.58 mm of the
# LA mid-surface, i.e. up to ~1.4 mm into the 4 mm LA); the softer penalty (L19_pen5) got past the control's wall. So the
# user's four-loft model (L18_lofts4_aggr: t = 1, 158 failed attempts, 1 h 36; its PVW_LA contact starts at t ~ 0.95)
# with the same softer penalty: is it faster? One change.
if 'L21_lofts4_pen' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L21_lofts4_pen')), 'L21_lofts4_pen exists; not overwriting'
    m = Model6(os.path.join(RUNS, 'L18_lofts4_aggr', 'L18_lofts4_aggr.feb'))
    assert m.set_contact('PVW_LA', penalty=0.5) == 1
    m.log.append('PVW_LA: penalty 5 -> 0.5 (auto_penalty 1 kept: 10x softer); the source contact pair is kept')
    emit('L21_lofts4_pen', m)

# L21D_allconn_pen, L21D_lofts4_pen (2026-09-25 ~03:30): the like-for-like runs L18D_allconn_aggr (t 0.985) and
# L18D_lofts4_aggr (t 0.987) crawled where 1-2 PVW_LA facets touch, as on the fast base; L19_pen5 (the softer penalty, the
# source contact kept) reached t = 1 there in 1 h 10 with 30 failed attempts. One change each: PVW_LA penalty 5 -> 0.5.
for _name, _base in (('L21D_allconn_pen', 'L18D_allconn_aggr'), ('L21D_lofts4_pen', 'L18D_lofts4_aggr')):
    if _name in WANT:
        assert not os.path.exists(os.path.join(RUNS, _name)), f'{_name} exists; not overwriting'
        m = Model6(os.path.join(RUNS, _base, _base + '.feb'))
        assert m.set_contact('PVW_LA', penalty=0.5) == 1
        m.log.append('PVW_LA: penalty 5 -> 0.5 (auto_penalty 1 kept: 10x softer); the source contact pair is kept')
        emit(_name, m)

# L22_allconn_pen_hold (2026-09-25 ~04:00): the equilibrium at the source loads on the model that passes the contact
# onset (L19_pen5: every loft as its connectors, real cutback, PVW_LA penalty 0.5; t = 1 in 1 h 10); L20_lofts4_hold
# (penalty 5) was crawling at the onset (t 0.989, 191 failed attempts). Held to t = 2.5, settle C 20 -> 2/s after t = 1
# (NOT IN SOURCE, after the ramp only).
if 'L22_allconn_pen_hold' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L22_allconn_pen_hold')), 'L22_allconn_pen_hold exists; not overwriting'
    m = Model6(os.path.join(RUNS, 'L19_pen5', 'L19_pen5.feb'))
    m.set_end_time(2.5)
    m.log.append('the load is held from t = 1 to 2.5 (Amp-1 extends CONSTANT)')
    bl = next(b for b in m.root.find('Loads') if b.get('name') == 'settle_damping')
    lc = next(l for l in m.root.find('LoadData') if l.get('id') == bl.find('C').get('lc'))
    assert [p.text for p in lc.find('points')] == ['1.0,0', '1.05,1'], [p.text for p in lc.find('points')]
    old = bl.find('C').text
    bl.find('C').text = '2'
    m.log.append(f'NOT IN SOURCE (after the ramp only): settle_damping C {old} -> 2/s, on at t = 1.0-1.05 (settle_on)')
    emit('L22_allconn_pen_hold', m)

# L23_lofts4_rho (2026-09-25 ~04:35): L21_lofts4_pen_r2 (the new best: t = 1, 9 failed attempts, 43 min) with its four
# lofts (AVW-Para L/R, CL L/R, USL L/R, PM: 1.29 kg at the old pipe-beam density 7.8e-7; the connectors they replace are
# massless) at tissue density 1.06e-9. One change: the faithful mass for the like-for-like dynamics. With all the lofts
# light, including the soft ones, the run was very slow (L15_rho: t 0.296 in 57 min).
if 'L23_lofts4_rho' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L23_lofts4_rho')), 'L23_lofts4_rho exists; not overwriting'
    m = Model6(os.path.join(RUNS, 'L21_lofts4_pen_r2', 'L21_lofts4_pen_r2.feb'))
    m.set_loft_density(1.06e-9)
    emit('L23_lofts4_rho', m)
