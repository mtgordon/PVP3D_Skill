"""Batch 43 (2026-09-25 night, user away ~8 h: "keep working on getting a bigger prolapse (especially posterior)").

Batch 42: the PVW at 50 % stiffness crawls at t ~ 0.66 whatever the apex support (CL/USL 50 % or 100 %, the USL-L loft on
shared nodes); there the PVW bulges forward (mean u_y +9.6, p90 +20.8 mm; the cervix +17.5) into the AVW and the canal's
inner contact (SlidingElastic1 = the source's Int-AVW-PVW-PEB-inner) is engaged over 580/884 and 1313/1508 facets. One
change each on L33_pvw50_clusl50 (lofts line, CL/USL 50 %, Vagina_PVW x 0.5, Load-LA 0.01, the rest 0.014), NOT IN SOURCE:
  L34_canalpen05        SlidingElastic1 penalty 5 -> 0.5 (auto_penalty kept; the source contact pair stays), as the PVW-LA
                        fix that got the source-load runs through
  L34_ramp3             the load ramped over 3 s (NOT IN SOURCE: the source ramps over 1 s)
  L34_canalpen05_ramp3  both (one change from each of the two above)
  L34_pvw70             the PVW at 70 % instead of 50 % (a milder softening)
usage: py -3.10 build_batch43.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L33_pvw50_clusl50', 'L33_pvw50_clusl50.feb')
WANT = set(sys.argv[1:])


def canal_pen(m):
    assert m.set_contact('SlidingElastic1', penalty=0.5) == 1
    m.log.append('SlidingElastic1 (the canal inner faces, AVW+cervix vs PVW+PeB): penalty 5 -> 0.5 (auto_penalty kept)')


def ramp3(m, s=3.0):
    ld = m.root.find('LoadData')
    for name in ('Amp-1', 'settle_on'):
        lc = next(l for l in ld if l.get('name') == name)
        pts = lc.find('points')
        old = [p.text for p in pts]
        for p in pts:
            t, v = p.text.split(',')
            p.text = f'{float(t) * s:g},{v}'
        m.log.append(f'load curve {name}: time x{s:g}: {old} -> {[p.text for p in pts]}')
    ctrl = m.root.find('Control')
    for path in ('step_size', 'time_stepper/dtmax'):
        el = ctrl.find(path)
        old = el.text
        el.text = '%g' % (float(old) * s)
        m.log.append(f'{path}: {old} -> {el.text}')
    m.log.append(f'NOT IN SOURCE: the load ramps over {s:g} s (Abaqus: 1 s); ends at t = '
                 f'{int(ctrl.find("time_steps").text) * float(ctrl.find("step_size").text):g}')


def pvw70(m):
    mat = next(x for x in m.root.find('Material') if x.get('name') == 'Vagina_PVW')
    for k in ('k', 'c1', 'c2'):
        mat.find(k).text = '%.6g' % (float(mat.find(k).text) * 0.7 / 0.5)
    m.log.append('NOT IN SOURCE (material): Vagina_PVW at 70 % of the original (was 50 %): c1, c2, k x 1.4')


VARIANTS = {'L34_canalpen05': (canal_pen,), 'L34_ramp3': (ramp3,), 'L34_canalpen05_ramp3': (canal_pen, ramp3),
            'L34_pvw70': (pvw70,)}
for name, fns in VARIANTS.items():
    if WANT and name not in WANT:
        continue
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(BASE)
    for fn in fns:
        fn(m)
    emit(name, m)

# L35_ramp3_pvw30 (2026-09-26 ~00:10): the 3 s ramp carries the soft PVW (L34_ramp3: 66 % load at 100 min, 6 failed attempts)
# where the 1 s runs crawl; for a bigger posterior prolapse, the PVW at 30 % of the original instead of 50 %. One change
# from L34_ramp3. NOT IN SOURCE.
if 'L35_ramp3_pvw30' in WANT:
    assert not os.path.exists(os.path.join(RUNS, 'L35_ramp3_pvw30')), 'L35_ramp3_pvw30 exists; not overwriting'
    m = Model6(os.path.join(RUNS, 'L34_ramp3', 'L34_ramp3.feb'))
    mat = next(x for x in m.root.find('Material') if x.get('name') == 'Vagina_PVW')
    for k in ('k', 'c1', 'c2'):
        mat.find(k).text = '%.6g' % (float(mat.find(k).text) * 0.3 / 0.5)
    m.log.append('NOT IN SOURCE (material): Vagina_PVW at 30 % of the original (was 50 %): c1, c2, k x 0.6')
    emit('L35_ramp3_pvw30', m)
