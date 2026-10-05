"""Batch 96 (2026-09-28 morning): seam springs instead of the seam tie (NOT IN SOURCE; the user's request).

The user (after the seam-tie results): "set up the seam springs runs. I don't like that the ties are affecting the results
that much." The tie holds the canal's lateral seam rigidly (116 AVW / cervix edge nodes to the PVW point each faces at
rest): on the Ogden-wall line it lifts the anterior wall 14.4 mm (L89_springs_newline_rhoi0_seamtie: Ba -17.7 vs -3.3).
Untied, that line's seam slides up to 16 mm and opens up to 21 mm; with the faithful walls up to 35 / 39 mm by t 0.76, and
those runs crawl. The springs join the same node / point pairs by a zero-length spring of k N/mm in every direction
(penalty-only linear constraints, build_batch90.seam_tie(maxaug=0)), so the edges can slide and open elastically.
k per node spans three decades around a rough estimate of the wall tissue's own shear stiffness across a ~3 mm fold,
G t / w x the ~1.6 mm node spacing: 0.02 (Yeoh walls, G ~0.012 MPa) to 0.16 N/mm (Ogden walls, G ~0.1 MPa).
One change each:
  L91_springs_newline_rhoi0_seamspr001 / _seamspr01 / _seamspr1   L87_springs_newline_rhoi0 (Ogden walls, untied) + springs
                                                                  of 0.01 / 0.1 / 1 N/mm: how the results move between
                                                                  free (Ba -3.3) and tied (-17.7)
  L91_newline_vwyeoh_rhoi0_seamspr001 / _seamspr01 / _seamspr1    L85_newline_vwyeoh_rhoi0 (the faithful walls, untied: crawled
                                                                  at t 0.78) + the same springs: do the faithful walls reach
                                                                  t = 1 with a softer seam?
usage: py -3.10 build_batch96.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit
from build_batch90 import seam_tie

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def springs(k):
    def fn(m):
        return seam_tie(m, penalty=k, maxaug=0)
    return fn


BUILDS = tuple((f'{pre}_seamspr{tag}', base, springs(k))
               for pre, base in (('L91_springs_newline_rhoi0', 'L87_springs_newline_rhoi0'),
                                 ('L91_newline_vwyeoh_rhoi0', 'L85_newline_vwyeoh_rhoi0'))
               for tag, k in (('001', 0.01), ('01', 0.1), ('1', 1)))
# 08:46: on the Ogden-wall line even 0.01 N/mm gives most of the tie's effect (Ba -14.6; free -3.3, tied -17.7; the seam
# slides 3.8 mm max instead of 16.4), so two softer ones find where the curve bends:
BUILDS += tuple((f'L91_springs_newline_rhoi0_seamspr{tag}', 'L87_springs_newline_rhoi0', springs(k))
                for tag, k in (('0001', 0.001), ('00001', 0.0001)))
# 09:20: on the Ogden-wall line 0.001 N/mm gives Ba -7.2 and 0.0001 gives -3.8 (free -3.3): do the faithful walls converge
# with springs that soft?
BUILDS += tuple((f'L91_newline_vwyeoh_rhoi0_seamspr{tag}', 'L85_newline_vwyeoh_rhoi0', springs(k))
                for tag, k in (('0001', 0.001), ('00001', 0.0001)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        print(name, ':', m.log[-1][:110])
        emit(name, m)
