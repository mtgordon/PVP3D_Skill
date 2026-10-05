"""Batch 36 (2026-09-25, user): add the Abaqus reference-only part PM_Plane to the FEBio models.

In PVP3DModel_job.inp, part PM_Plane (396 nodes, 350 S4R) is instanced as PM_Plane-1 with no transform and marked
*Display Body (Constraint-PM_Plane): shown in Abaqus but not part of the analysis (no section, no material, no tie or
connector). FEBio has no display body, so it is added as a separate FEBio rigid body fixed in all 6 DOFs: its own new
nodes (shared with nothing), quad4 elements, a "rigid body" material and a rigid_fixed constraint. It adds no unknowns
and cannot change the solution; it only shows in FEBio Studio. Shell thickness 0.5 mm for display only (the source part
has no section). Syntax verified in claude_diag/mini/rigid_display.feb (rigid shell nodes stay at u = 0; a rigid shell
domain takes no shell_normal_nodal).
  L26_lofts_la3_pm    L25_lofts_la3_opt10 (the lofts line, current fastest) + PM_Plane
  L26_springs_la3_pm  L25_springs_la3_opt10 (the springs line = the Abaqus replica line, current fastest) + PM_Plane
  L26D_allconn_pen_pm L21D_allconn_pen (the Abaqus U/V comparison model) + PM_Plane (built, not run)
The user wants PM_Plane in every new model from now on (add_display_body).
usage: py -3.10 build_batch36.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants6 import Model6, JOBS, emit
from paths import INP_FILE

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASES = {'L26_lofts_la3_pm': 'L25_lofts_la3_opt10', 'L26_springs_la3_pm': 'L25_springs_la3_opt10',
         'L26D_allconn_pen_pm': 'L21D_allconn_pen'}   # the Abaqus-comparison model + PM_Plane (user: build, not run)
WANT = set(sys.argv[1:])


def read_part(inp, part):
    """Nodes and S4R elements of one Abaqus part, and a check that its instance has no transform."""
    lines = open(inp).read().splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.strip().lower() == f'*part, name={part.lower()}')
    nodes, elems, mode = {}, [], None
    for l in lines[i0 + 1:]:
        s = l.strip()
        if s.lower().startswith('*end part'):
            break
        if s.startswith('*'):
            mode = 'n' if s.lower() == '*node' else ('e' if s.lower().startswith('*element, type=s4r') else None)
            assert mode or not s.lower().startswith('*element'), s
            continue
        v = [x.strip() for x in s.split(',') if x.strip()]
        if mode == 'n':
            nodes[int(v[0])] = tuple(float(x) for x in v[1:4])
        elif mode == 'e':
            elems.append(tuple(int(x) for x in v[1:5]))
    j = next(i for i, l in enumerate(lines) if l.lower().startswith(f'*instance, name={part.lower()}-1'))
    assert lines[j + 1].strip().lower() == '*end instance', 'the instance has a transform: apply it first'
    return nodes, elems


def add_display_body(m, part='PM_Plane', thickness=0.5, rho=1.06e-9):
    nodes, elems = read_part(INP_FILE, part)
    order = sorted(nodes)
    new_ids = m.add_nodes(part, [nodes[a] for a in order])
    fid = dict(zip(order, new_ids))
    mesh = m.root.find('Mesh')
    blocks = mesh.findall('Elements')
    eid = max(int(e.get('id')) for b in blocks for e in b) + 1
    blk = ET.Element('Elements', {'type': 'quad4', 'name': part})
    for k, e in enumerate(elems):
        ET.SubElement(blk, 'elem', {'id': str(eid + k)}).text = ','.join(str(fid[a]) for a in e)
    mesh.insert(list(mesh).index(blocks[-1]) + 1, blk)
    mats = m.root.find('Material')
    mid = max(int(x.get('id')) for x in mats) + 1
    mname = f'{part}_display'
    mat = ET.SubElement(mats, 'material', {'id': str(mid), 'name': mname, 'type': 'rigid body'})
    ET.SubElement(mat, 'density').text = '%g' % rho
    dom = ET.SubElement(m.root.find('MeshDomains'), 'ShellDomain', {'name': part, 'mat': mname})
    ET.SubElement(dom, 'shell_thickness').text = '%g' % thickness
    rigid = m.root.find('Rigid')
    if rigid is None:
        kids = list(m.root)
        rigid = ET.Element('Rigid')
        m.root.insert(kids.index(m.root.find('Boundary')) + 1, rigid)
    bc = ET.SubElement(rigid, 'rigid_bc', {'name': f'{part}_fixed', 'type': 'rigid_fixed'})
    ET.SubElement(bc, 'rb').text = mname
    for d in ('Rx', 'Ry', 'Rz', 'Ru', 'Rv', 'Rw'):
        ET.SubElement(bc, f'{d}_dof').text = '1'
    m.log.append(f'added the Abaqus display body {part} (*Display Body in the source: reference only, no section, no '
                 f'material, no connections): {len(nodes)} nodes (new ids {new_ids[0]}-{new_ids[-1]}, shared with '
                 f'nothing), {len(elems)} S4R -> quad4 (ids {eid}-{eid + len(elems) - 1}), rigid body material {mname} '
                 f'(density {rho:g}) fixed in all 6 DOFs (rigid_fixed): adds no unknowns, cannot change the solution. '
                 f'Shell thickness {thickness:g} mm for display only (the source part has none)')


for name, base in BASES.items():
    if WANT and name not in WANT:
        continue
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    m = Model6(os.path.join(RUNS, base, base + '.feb'))
    add_display_body(m)
    emit(name, m)
