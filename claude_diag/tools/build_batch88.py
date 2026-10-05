"""Batch 88 (2026-09-27 ~20:10): keep pushing the faithful walls (the user, ~20:00) on the new springs line.

L84_springs_newline_vwyeoh = the new springs line (the 30 pinch facets out, PVW_LA penalty 0.5 + seg_up 2, the healthy LA,
Load-LA 1/3, the PeB refit) + the Yeoh vaginal walls. With PVW_LA penalty 0.5 the walls crawled past t ~0.65 at the canal
seam (L79_springs_la3_vwyeoh_pvwpen05). Canal settings never tried together with the PVW_LA fix, one change each:
  L85_newline_vwyeoh_se1pen05     SlidingElastic1 penalty 5 -> 0.5
  L85_newline_vwyeoh_rhoi0        solver: rhoi 0.5 -> 0 (full high-frequency damping; the seam nodes reverse each step)
  L85_newline_vwyeoh_se1stol01    SlidingElastic1 search_tol 0.01 -> 0.1 (a node just off the last facet still counts)
  L85_newline_vwyeoh_se1knmult1   SlidingElastic1 knmult 0 -> 1 (the exact contact stiffness; the forum's advice)
FEBio 4.13's "contact potential" was tried first on a mini model (claude_diag/contact_potential/cp_mini.py): it failed as
soon as the blocks slid (kc 3e-5 .. 1e-2, R_in 0.05, R_out 0.1), while sliding-elastic ran cleanly; parked.
usage: py -3.10 build_batch88.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch79 import se1, rhoi0

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])

BUILDS = (('L85_newline_vwyeoh_se1pen05', se1(penalty=0.5)),
          ('L85_newline_vwyeoh_rhoi0', rhoi0),
          ('L85_newline_vwyeoh_se1stol01', se1(search_tol=0.1)),
          ('L85_newline_vwyeoh_se1knmult1', se1(knmult=1)))
if __name__ == '__main__':
    base = 'L84_springs_newline_vwyeoh'
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
