"""Batch 76 (2026-09-27): the LA on FEBio's old (pre-2.6) mid-surface shell formulation, the closest FEBio has to the
Abaqus S4R (nodes on the mid-surface, directors free at pinned nodes).

The user's question (2026-09-27): is the LA thickness treated the same in both codes? FEBio 4.13 still has the old
formulation as the ShellDomain type "elastic-shell-old" (found in febiomech.dll). It rejects <shell_thickness> in the
domain; the thickness is given per element in <MeshData><ElementData type="shell thickness"> (verified on the panel test
claude_diag/la_stiffness/panel2_yeoh_shell_old_md.feb: 30.57 mm centre deflection against 30.67 for elastic-shell).
It takes no shell_normal_nodal (averaged nodal normals, extending +-t/2).
  L69_pen5_laold   L19_pen5 with the 4 LA ShellDomains as elastic-shell-old, thickness 4.0 per element. NOT IN SOURCE
                   (a formulation change: the source is an S4R, which this approximates)
usage: py -3.10 build_batch76.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def laold(m):
    doms = m.root.find('MeshDomains')
    eb = m.elem_blocks()
    md = m.root.find('MeshData')
    if md is None:
        md = ET.Element('MeshData')
        kids = list(m.root)
        m.root.insert(kids.index(doms) + 1, md)
    done = []
    for d in doms:
        name = d.get('name', '')
        if not name.startswith('LA_'):
            continue
        t = d.find('shell_thickness').text
        for tag in ('shell_thickness', 'shell_normal_nodal'):
            el = d.find(tag)
            if el is not None:
                d.remove(el)
        d.set('type', 'elastic-shell-old')
        blk = eb[name]
        nn = 3 if blk.get('type') == 'tri3' else 4
        ed = ET.SubElement(md, 'ElementData', {'name': name + '_thickness', 'type': 'shell thickness', 'elem_set': name})
        for i in range(1, len(blk) + 1):
            ET.SubElement(ed, 'e', {'lid': str(i)}).text = ','.join([t] * nn)
        done.append(f'{name} ({len(blk)} {blk.get("type")}, {t} mm)')
    assert len(done) == 4, done
    m.log.append('NOT IN SOURCE (formulation; the source LA is an Abaqus S4R): the LA on FEBio\'s old mid-surface shell '
                 'formulation, ShellDomain type elastic-shell-old, thickness per element in MeshData (shell_thickness and '
                 'shell_normal_nodal removed from the domains): ' + ', '.join(done))


BUILDS = (('L69_pen5_laold', 'L19_pen5', laold),)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
