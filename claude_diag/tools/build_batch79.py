"""Batch 79 (2026-09-27 ~08:55): the canal-contact crawl with the source-faithful (Yeoh) vaginal walls.

Diagnosis (tools/contact_onset.py on L73_pen5_paperP1_vwyeoh_nopinch t 0.555-0.585 and L70_pen5_vwyeoh t 0.615-0.63):
the fastest nodes are AVW nodes on the edge row of SlidingElastic1Secondary (2069-2072, 2079, 2159-2162), reversing
direction between converged steps at 0.1-2 m/s, with pressure only from the two-pass second pass (AVW points projected
onto PVW facets). The AVW and PVW share no nodes: the canal is closed laterally by this contact alone, and 145 of the
secondary's 168 edge nodes lie within 0.75 mm of the primary's edge (a seam along the whole length). A node there drops on
and off the PVW's last facet. FEBio manual 3.13.1: facet switching can leave "a node [that] oscillates continuously
between two adjacent facets and thus prevents FEBio from meeting the displacement convergence tolerance" (the failed
attempts meet the energy norm, never the displacement norm); seg_up limits the updates.
One change each from L74_nopinch_vwyeoh (= L71_aggr_noparcusfacets + the Yeoh walls, source loads):
  L75_nopinch_vwyeoh_se1segup2    SlidingElastic1 seg_up 0 -> 2 (facet switching only in each step's first 2 iterations)
  L75_nopinch_vwyeoh_se1pen05     SlidingElastic1 penalty 5 -> 0.5 (auto_penalty kept: 10x softer; as PVW_LA's Test B fix)
  L75_nopinch_vwyeoh_se1twopass0  SlidingElastic1 two_pass 1 -> 0 (the rattle rides the second pass)
  L75_nopinch_vwyeoh_rhoi0        solver only: rhoi 0.5 -> 0 (the generalized-alpha integrator's full high-frequency
                                  damping; the edge nodes reverse between converged steps, a period-2 oscillation)
  L75_nopinch_vwyeoh_seamsec      NOT IN SOURCE: the 142 SlidingElastic1Secondary facets (AVW / cervix) holding a node
                                  within 0.75 mm of the primary's edge removed (the seam), so the AVW's contact edge sits
                                  ~one row inside the PVW's surface instead of on its edge
usage: py -3.10 build_batch79.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import collections
import os
import sys

import numpy as np

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def se1(**kw):
    def fn(m):
        assert m.set_contact('SlidingElastic1', **kw) == 1
    return fn


def rhoi0(m):
    m.set_control(rhoi=0)
    m.log.append('solver only: rhoi 0.5 -> 0 (generalized-alpha, full high-frequency damping); the model is unchanged')


def _bedges(fs):
    c = collections.Counter()
    for f in fs:
        for i in range(len(f)):
            c[tuple(sorted((f[i], f[(i + 1) % len(f)])))] += 1
    return [e for e, k in c.items() if k == 1]


def _dseg(p, a, b):
    ab = b - a
    t = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0, 1)
    return np.linalg.norm(p - (a + t * ab))


def seam_out(m, side='Secondary', other='Primary', tol=0.75, pair='SlidingElastic1'):
    """NOT IN SOURCE: drop the <pair><side> facets holding a node within tol of the <pair><other> surface's edge."""
    X = m.nodes()
    surf = {s.get('name'): s for s in m.mesh.findall('Surface')}
    S, O = surf[pair + side], surf[pair + other]
    fl = lambda s: [tuple(int(v) for v in f.text.split(',')) for f in s]
    eS, eO = _bedges(fl(S)), _bedges(fl(O))
    bS = {v for e in eS for v in e}
    seam = {v for v in bS if min(_dseg(X[v], X[a], X[b]) for a, b in eO) < tol}
    n0 = len(S)
    gone = [f for f in list(S) if {int(v) for v in f.text.split(',')} & seam]
    for f in gone:
        S.remove(f)
    for i, f in enumerate(S, 1):
        f.set('id', str(i))
    m.log.append(f'NOT IN SOURCE (canal-contact diagnostic): {len(gone)} of {n0} {pair}{side} facets removed: those holding '
                 f'one of the {len(seam)} edge nodes within {tol} mm of the {pair}{other} edge (the lateral seam where the '
                 f'AVW and PVW contact surfaces meet edge to edge); {len(S)} left')
    return len(gone)


BUILDS = (('L75_nopinch_vwyeoh_se1segup2', se1(seg_up=2)),
          ('L75_nopinch_vwyeoh_se1pen05', se1(penalty=0.5)),
          ('L75_nopinch_vwyeoh_se1twopass0', se1(two_pass=0)),
          ('L75_nopinch_vwyeoh_rhoi0', rhoi0),
          ('L75_nopinch_vwyeoh_seamsec', seam_out))
if __name__ == '__main__':
    base = 'L74_nopinch_vwyeoh'
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
