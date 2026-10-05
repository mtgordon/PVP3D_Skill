"""Batch 44 (2026-09-26, user away ~8 h): what keeps the perineal body (_PickedSet66) from rotating about x.

User: "I expect PickedSet66 to rotate a lot more ... see what is preventing it from rotating about the x axis. Try not to
weaken the vaginal wall material (AVW, PVW, Cervix, and even PeB). Try adjusting the connective tissue instead. Ideally
this would work for both the lofted and the non-lofted versions." Wanted sense: the cranial end +y and -z, up to 90 deg
(tools/peb_rotation.py rot_-x > 0). Bases: L26_springs_la3_pm and L26_lofts_la3_pm (user). Diagnostic runs may release
anything, one thing at a time, labelled DIAGNOSTIC; the proposed fixes change connective tissue only (user).

In the bases the body turns the wanted way up to t ~0.66 (lofts +4.6, springs +9.4 deg), then turns back (t = 1: -2.5 /
+3.4 deg) while its centroid stops descending. Its attachments (tools/what_holds.py): 108 nodes shared with the PVW; the
sphincter connectors to LA_PCMPRM (16 side, 24 posterior); PM_PeB (8 + 8) and PeB-constrin (26) to fixed points, which
carry ~0 (impaired, as in the source); 2 of the 26 P-arcus connectors; BC-VW-mid (x only, no moment about x); its top face
in the canal contact (not touching) and under Load-PeB-top.

DIAGNOSTIC cuts, one each, on both lines (a spring set cut = its material scale 0, so the dmat list positions stay):
  _nosphside   LA_sphincter_side_conn (16) cut
  _nosphpost   LA_sphincter_post_conn (24) cut
  _noparcus    Parcus_conn (26, the PVW's posterior-arcus connectors) cut
  _nopvwla     the PVW-LA contact removed (the PVW's back face against the LA)
usage: py -3.10 build_batch44.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASES = {'springs': os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb'),
         'lofts': os.path.join(RUNS, 'L26_lofts_la3_pm', 'L26_lofts_la3_pm.feb')}
WANT = set(sys.argv[1:])


def spring_scale(m, set_name, s, label='DIAGNOSTIC'):
    disc = m.root.find('Discrete')
    d = next(e for e in disc.findall('discrete') if e.get('discrete_set') == set_name)
    dm = disc.findall('discrete_material')[int(d.get('dmat')) - 1]
    el = dm.find('scale')
    old = el.text
    el.text = '%g' % (float(old) * s)
    n = sum(1 for ds in m.root.find('Mesh').findall('DiscreteSet') if ds.get('name') == set_name
            for _ in ds.findall('delem'))
    m.log.append(f'{label}: {set_name} ({n} springs), material {dm.get("name")} scale {old} -> {el.text}')


def nosphside(m):
    spring_scale(m, 'LA_sphincter_side_conn', 0.0)


def nosphpost(m):
    spring_scale(m, 'LA_sphincter_post_conn', 0.0)


def noparcus(m):
    spring_scale(m, 'Parcus_conn', 0.0)


def nopvwla(m):
    m.remove_contact('PVW_LA')
    m.log.append('DIAGNOSTIC: the PVW-LA contact (the source contact pair) removed')


VARIANTS = {'nosphside': nosphside, 'nosphpost': nosphpost, 'noparcus': noparcus, 'nopvwla': nopvwla}
if __name__ == '__main__':
    for line, base in BASES.items():
        for tag, fn in VARIANTS.items():
            name = f'L36_{line}_{tag}'
            if WANT and name not in WANT:
                continue
            assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
            m = Model6(base)
            fn(m)
            emit(name, m)
