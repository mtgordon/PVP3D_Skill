"""Batch 78 (2026-09-27 ~08:45): the user's PVW_LA choice on both lines, the faithful materials with it, and the first
ballooning screens.

The user (2026-09-27 08:40): "For now, let's remove the 30 facets. Use the spring lines for comparison, but you can use
the lofts for faster runs." = option (b) of NEXT_SESSION_PROMPT_2026-09-27.md: the 30 PVW_LA primary facets that hold a
P-arcus connector end (a Parcus_conn node) are removed from the contact surface and the source-scale penalty 5 is kept
(L71_aggr_noparcusfacets: t = 1 with the same solution as penalty 0.5 within 0.1 mm). NOT IN SOURCE.

  L74_springs_la3_nopinch  L26_springs_la3_pm (the springs line) + option (b): PVW_LA penalty 0.5 -> 5, 30 facets out
  L74_lofts_la3_nopinch    L26_lofts_la3_pm (the lofts line) + the same
Task 1, the source's tissue laws refit (Yeoh, as the LA already is), at the source loads, one change each from
L71_aggr_noparcusfacets (= L19_pen5, the Abaqus replica, with option (b)):
  L74_nopinch_vwyeoh       vaginal walls (AVW, PVW, cervix) Ogden -> Yeoh c1 0.005964127 c2 0.01883547 k 1
  L74_nopinch_pebyeoh      perineal body Ogden -> Yeoh c1 0.05106451 c2 0.1059592 k 5.7
Task 2, the LA's ballooning (area x2.16 at the source loads), one change each from L71_aggr_noparcusfacets, NOT IN SOURCE:
  L74_nopinch_la100        the healthy LA: the source's own PCM-LA_Yamada100% table (defined, unused), fit_yeoh.py:
                           c1 0.02347128 c2 0.03283498, k 2.0 (nu 0.47 at ~40 % strain, the rule used for the 50 % law);
                           about 2x the 50 % law above 30 % strain, 2.4x at 10 %
  L74_nopinch_la0          no pressure on the LA (Load-LA 0.014 -> 0): how much of the stretch the organs alone cause
usage: py -3.10 build_batch78.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch77 import vwyeoh, pebyeoh

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
LA100 = {'c1': '0.02347128', 'c2': '0.03283498', 'k': '2'}


def nopinch(m, penalty='5'):
    """Option (b): drop every PVW_LA primary facet holding a P-arcus connector end; PVW_LA penalty back to 5."""
    mesh = m.mesh
    pn = {int(v) for d in mesh.findall('DiscreteSet') if d.get('name') == 'Parcus_conn' for e in d
          for v in e.text.split(',')}
    assert len(pn) > 0
    s = next(s for s in mesh.findall('Surface') if s.get('name') == 'PVW_LA_primary')
    n0 = len(s)
    gone = [f for f in list(s) if {int(v) for v in f.text.split(',')} & pn]
    for f in gone:
        s.remove(f)
    for i, f in enumerate(s, 1):
        f.set('id', str(i))
    assert (n0, len(gone)) == (676, 30), (n0, len(gone))
    c = next(c for c in m.root.find('Contact') if c.get('name') == 'PVW_LA')
    old = c.find('penalty').text
    c.find('penalty').text = penalty
    m.log.append(f'NOT IN SOURCE (the user\'s choice (b), 2026-09-27): the {len(gone)} PVW_LA primary facets holding a '
                 f'P-arcus connector end (Parcus_conn) removed from the contact surface ({len(s)} of {n0} left), as in '
                 f'L71_aggr_noparcusfacets; PVW_LA penalty {old} -> {penalty} (auto_penalty 1)')


def la100(m):
    mat = m._material('LA_Yamada50pct_Yeoh')
    old = {k: mat.find(k).text for k in ('c1', 'c2', 'k')}
    for k, v in LA100.items():
        mat.find(k).text = v
    m.log.append(f'NOT IN SOURCE (ballooning screen): the healthy LA, the source\'s own PCM-LA_Yamada100% table (defined in '
                 f'the .inp, unused) fit to Yeoh (fit_yeoh.py; k = nu 0.47 at ~40 % strain): LA_Yamada50pct_Yeoh {old} -> '
                 f'{LA100} (material name kept)')


def la0(m):
    m.scale_surface_load('Load-LA', 0.0)


BUILDS = (('L74_springs_la3_nopinch', 'L26_springs_la3_pm', nopinch),
          ('L74_lofts_la3_nopinch', 'L26_lofts_la3_pm', nopinch),
          ('L74_nopinch_vwyeoh', 'L71_aggr_noparcusfacets', vwyeoh),
          ('L74_nopinch_pebyeoh', 'L71_aggr_noparcusfacets', pebyeoh),
          ('L74_nopinch_la100', 'L71_aggr_noparcusfacets', la100),
          ('L74_nopinch_la0', 'L71_aggr_noparcusfacets', la0))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
