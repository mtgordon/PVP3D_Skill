"""Batch 74 (2026-09-27, the user away ~8 h): what sets the LA displacement at the source loads, and two contact settings
for the PVW_LA stall.

LA question (the user, 2026-09-27: why does the LA move so much more in FEBio than in Abaqus; is the thickness treated
differently?). At t = 1 of L19_pen5 the LA's tied edge (the 199 LA nodes held to the ATLA / posterior-arcus truss chains)
moves as much as the rest (median 33.6 vs 33.1 mm), the chains at +17..52 % strain: the LA hangs on the chains
(tools/la_support_split.py). One change each from L19_pen5 (the Abaqus replica at the source loads):
  L67_pen5_ctrl      identical copy (control, same batch)
  L67_pen5_laback    the LA's pinned BCs hold the whole section: the back face (shell displacement sx, sy, sz) fixed
                     too on BC-LA-ICM-posterior / PCM-anterior / PRM-anterior, and sx on BC-LA-mid. FEBio pins only the
                     front face; an Abaqus S4R PINNED node holds the whole section's translation. Holding both faces
                     also stops the section rotating there (clamped), so this bounds the effect from above. NOT IN SOURCE
  L67_pen5_lat8      LA shell thickness 4 -> 8 mm (sensitivity: membrane x2, bending x8). NOT IN SOURCE
  L67_pen5_chain10   the ATLA and posterior-arcus truss-chain springs 10x stiffer (sensitivity). NOT IN SOURCE
PVW_LA stall (Test B): the LA touches the lateral edge of the PVW contact surface (the P-arcus connector ends), a
projection there drops on and off the edge row, and the edge nodes rattle. One change each from L18_allconn_aggr
(penalty 5, stalled at t 0.9036-0.914):
  L67_aggr_knmult1   PVW_LA knmult 0 -> 1 (the full contact stiffness; the forum: always recommended)
  L67_aggr_stol01    PVW_LA search_tol 0.01 -> 0.1 (projections just past a facet edge still count as inside)
usage: py -3.10 build_batch74.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
LA_PINS = ('BC-LA-ICM-posterior', 'BC-LA-PCM-anterior', 'BC-LA-PRM-anterior')


def control(m):
    m.log.append('identical copy (the control for batch 74)')


def laback(m):
    bnd = m.root.find('Boundary')
    for ns, dofs in [(n, (1, 1, 1)) for n in LA_PINS] + [('BC-LA-mid', (1, 0, 0))]:
        m._nodeset(ns)
        bc = ET.SubElement(bnd, 'bc', {'name': ns + '_back', 'node_set': ns, 'type': 'zero shell displacement'})
        for tag, v in zip(('sx_dof', 'sy_dof', 'sz_dof'), dofs):
            ET.SubElement(bc, tag).text = str(v)
        m.log.append(f'NOT IN SOURCE (FEBio shells pin only the front face; an Abaqus S4R PINNED node holds the whole '
                     f'section): {ns} back face held too, zero shell displacement sx/sy/sz = {dofs}')


def lat8(m):
    n = 0
    for d in m.root.find('MeshDomains'):
        if d.get('name', '').startswith('LA_'):
            d.find('shell_thickness').text = '8.0'
            n += 1
    assert n == 4, n
    m.log.append('NOT IN SOURCE (sensitivity): LA shell thickness 4.0 -> 8.0 in the 4 LA ShellDomains')


def chain10(m):
    names = ('ATLA-Hyper_truss_mat', 'Posterior_Arcus-Hyper_truss_mat')
    for mat in m.root.find('Discrete').findall('discrete_material'):
        if mat.get('name') in names:
            mat.find('scale').text = '10'
    m.log.append(f'NOT IN SOURCE (sensitivity): {names} scale 1 -> 10 (the truss chains 10x stiffer)')


def contact(**kw):
    def fn(m):
        assert m.set_contact('PVW_LA', **kw) == 1
    return fn


BUILDS = (('L67_pen5_ctrl', 'L19_pen5', control),
          ('L67_pen5_laback', 'L19_pen5', laback),
          ('L67_pen5_lat8', 'L19_pen5', lat8),
          ('L67_pen5_chain10', 'L19_pen5', chain10),
          ('L67_aggr_knmult1', 'L18_allconn_aggr', contact(knmult=1)),
          ('L67_aggr_stol01', 'L18_allconn_aggr', contact(search_tol=0.1)))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
