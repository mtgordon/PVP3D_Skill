"""Batch 80 (2026-09-27 ~09:05): more LA ballooning screens (task 2), one change each from L71_aggr_noparcusfacets
(the Abaqus replica with the user's option (b), source loads; LA area x2.15 at t = 1: ICM x2.29, PCM x1.83, PCMPRM x1.84).
  L76_nopinch_chainpin  NOT IN SOURCE: every ATLA / posterior-arcus chain node held (zero displacement), as if the arcus
                        were attached along the pelvic sidewall; the source pins each chain at its two ends only and the
                        chains start 26 % longer than their chords. Separates the LA's own stretch from the supports'
                        sag (the tied edge stretches 1.59 at the source loads)
  L76_nopinch_la200     NOT IN SOURCE: the LA at 2x the healthy law (PCM-LA_Yamada100% x 2, Yeoh c1 0.04694256
                        c2 0.06566996 k 4): with L74_nopinch_la100 and the 50 % law, how stretch falls with stiffness
usage: py -3.10 build_batch80.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])
CHAINS = ('ATLA_Left_springs', 'ATLA_Right_springs', 'Posterior_Arcus_Left_springs', 'Posterior_Arcus_Right_springs')


def chainpin(m):
    nodes = sorted({int(v) for d in m.mesh.findall('DiscreteSet') if d.get('name') in CHAINS for e in d
                    for v in e.text.split(',')})
    assert len(nodes) == 11 * 2 + 18 * 2 + 4, len(nodes)   # 4 chains: springs + 1 nodes each
    ns = ET.SubElement(m.mesh, 'NodeSet', {'name': 'chain_all_pinned'})
    ns.text = ', '.join(str(n) for n in nodes)
    bnd = m.root.find('Boundary')
    bc = ET.SubElement(bnd, 'bc', {'name': 'chain_all_pinned', 'node_set': 'chain_all_pinned', 'type': 'zero displacement'})
    for d in ('x_dof', 'y_dof', 'z_dof'):
        ET.SubElement(bc, d).text = '1'
    m.log.append(f'NOT IN SOURCE (ballooning screen): all {len(nodes)} nodes of the 4 arcus chains ({", ".join(CHAINS)}) '
                 f'held in x, y, z (the source pins each chain at its 2 ends only)')


def la200(m):
    mat = m._material('LA_Yamada50pct_Yeoh')
    old = {k: mat.find(k).text for k in ('c1', 'c2', 'k')}
    new = {'c1': '0.04694256', 'c2': '0.06566996', 'k': '4'}
    for k, v in new.items():
        mat.find(k).text = v
    m.log.append(f'NOT IN SOURCE (ballooning screen): the LA at 2x the healthy law (the source\'s PCM-LA_Yamada100% Yeoh fit '
                 f'x 2): LA_Yamada50pct_Yeoh {old} -> {new} (material name kept)')


BUILDS = (('L76_nopinch_chainpin', chainpin), ('L76_nopinch_la200', la200))
if __name__ == '__main__':
    base = 'L71_aggr_noparcusfacets'
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        emit(name, m)
