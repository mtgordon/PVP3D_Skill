"""Batch 77 (2026-09-27 ~03:00): the tissue materials refit to the source's Marlow laws, and Test B's next contact setting.

Found 03:00: the FEBio vaginal-wall and perineal-body materials are Ogden fits that do not reproduce the source's Marlow
*Uniaxial Test Data (skill scripts/fit_yeoh.py, --ogden comparison, nominal stress, incompressible):
  Vagina_AVW  Ogden c1 0.1 m1 1, c2 0.02 m2 5.159, k 100:     uniaxial 3x stiffer at 20 % strain; equibiaxial 2.1x softer
                                                               at stretch 1.5 (0.091 vs 0.191 MPa)
  Vagina_PVW / Vagina_Cervix  Ogden c1 0.0887 m1 1, c2 0.0193 m2 5.159, k 100: equibiaxial 2.3x softer at 1.5
  PeB-Vagina500%stiffer  Ogden c1 0.0641 m1 5.97, c2 0.0318 m2 -8.04, k 100: uniaxial 2.9-4.1x softer (0.092 vs 0.38 MPa
                                                               at 1.5), planar 2.8-3.2x softer
The LA had the same trap and was refit to Yeoh on 2026-09-23 (gotcha 23: Marlow is I1-only, so is Yeoh). The Yeoh fits
(fit_yeoh.py) match Marlow in uniaxial, planar and equibiaxial stretch. k: the value that reproduces the source's
constant-nu (0.47) volume change at ~40 % uniaxial strain, as chosen for the LA (vagina 0.943 -> 1.0; PeB 5.66 -> 5.7);
the current walls use k 100. One change each (these ARE the source's laws, refit; not a model change):
  L70_pen5_vwyeoh            L19_pen5, Vagina_AVW, Vagina_PVW, Vagina_Cervix -> Yeoh c1 0.005964127 c2 0.01883547 k 1
  L70_pen5_pebyeoh           L19_pen5, PeB-Vagina500%stiffer -> Yeoh c1 0.05106451 c2 0.1059592 k 5.7
  L70_pen5_paperP1_vwyeoh    L68_pen5_paperP1 (the paper's case P1) + the vaginal walls as above
  L70_pen5_paperP2_vwyeoh    L68_pen5_paperP2 (the paper's case P2) + the vaginal walls as above (built inline 03:10)
Test B (the PVW_LA stall; knmult 1 and search_tol 0.1 did not help):
  L70_aggr_segup2            L18_allconn_aggr, PVW_LA seg_up 0 -> 2 (built inline 03:00; each contact point keeps its facet after 2 iterations)
  L70_aggr_twopass0          L18_allconn_aggr, PVW_LA two_pass 1 -> 0 (the second pass, LA nodes onto the PVW edge
                             facets, carried all of the contact at the onset; the first pass never touched)
usage: py -3.10 build_batch77.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
VW = {'c1': '0.005964127', 'c2': '0.01883547', 'k': '1'}
PEB = {'c1': '0.05106451', 'c2': '0.1059592', 'k': '5.7'}


def to_yeoh(m, name, p, why):
    mat = m._material(name)
    old = {c.tag: c.text for c in mat}
    dens = mat.find('density').text
    for c in list(mat):
        mat.remove(c)
    mat.set('type', 'Yeoh')
    ET.SubElement(mat, 'density').text = dens
    for k in ('c1', 'c2', 'k'):
        ET.SubElement(mat, k).text = p[k]
    keep = {k: old[k] for k in ('c1', 'm1', 'c2', 'm2', 'k') if k in old}
    m.log.append(f'the source\'s Marlow law refit ({why}): {name} Ogden {keep} -> Yeoh {p} (fit_yeoh.py on the source\'s '
                 f'*Uniaxial Test Data; k reproduces the source\'s nu 0.47 volume change at ~40 % strain)')


def vwyeoh(m):
    for n in ('Vagina_AVW', 'Vagina_PVW', 'Vagina_Cervix'):
        to_yeoh(m, n, VW, 'the Ogden fit was up to 2.3x softer in equibiaxial stretch and 3x stiffer at small strain')


def pebyeoh(m):
    to_yeoh(m, 'PeB-Vagina500%stiffer', PEB, 'the Ogden fit was 2.9-4.1x softer in uniaxial stretch')


def twopass0(m):
    assert m.set_contact('PVW_LA', two_pass=0) == 1


BUILDS = (('L70_pen5_vwyeoh', 'L19_pen5', vwyeoh),
          ('L70_pen5_paperP1_vwyeoh', 'L68_pen5_paperP1', vwyeoh),
          ('L70_pen5_pebyeoh', 'L19_pen5', pebyeoh),
          ('L70_aggr_twopass0', 'L18_allconn_aggr', twopass0))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
