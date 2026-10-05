"""Batch 30 (2026-09-25, overnight): what stalls the stab-spring-free model at the source pressures (batch 29).

Batch 29 (build_batch29.py, base "src" = L11_stab with Load-LA back at 0.014 MPa): every 1 s run crawled at t ~ 0.62-0.66
(~72 % load) whatever the one change (exact pressure tangent, AVW-Para or CL/USL lofts as connectors, lofts at tissue
density, mass damping, damping + 3 s ramp); only L15_conn (all 12 non-P-arcus lofts as their connectors) got further
(0.758). At the wall the fastest nodes (270-600 mm/s) and the lowest J were in the near-zero PM_PeB lofts (c1 4e-5 MPa at
tissue density: almost no stiffness and no mass), except in L15_conn, which has none; many retries then failed in their
first iteration (the line search cut the Newton update to 6 %, then 115 negative jacobians), and with aggressiveness 0
each retry takes only dt0/21 off the step. One change each:
  L16_nzconn       src + the five near-zero lofts (PM_PeB L/R, PM_avw_bottom L/R, PeB-constrin) -> their 46 Abaqus
                   connectors as springs (as L7_nzconn; in the source they carry ~0.003 N at 47 mm); the user's other
                   lofts stay
  L16_nzconn_aggr  L16_nzconn + a real cutback (aggressiveness 0 -> 1, cutback 0.25 -> 0.5: halve the step at each
                   failure). Solver only
  L16_aggr         src + the same cutback. Solver only
  L16_allconn      L15_conn + the two P-arcus lofts -> their 26 connectors (every Abaqus connector as a spring, no
                   lofts: the source's supports exactly; LPC_DM1_pconn had them with the stab springs)
  L16_allconn2     L16_allconn fixed: it failed its first step (850 negative jacobians at ~zero load) because the 8 + 8
                   loft-edge nodes that batch 19 put into the posterior-arcus spring chains have no element and no mass
                   once the P-arcus lofts are gone (singular); unsplit_chain() merges each one's two springs back
  L16_conn_aggr    L15_conn + the cutback. Solver only
usage: py -3.10 build_batch30.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys
import xml.etree.ElementTree as ET

from variants6 import Model6, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = os.path.join(RUNS, 'L11_stab', 'L11_stab.feb')
CONN = os.path.join(RUNS, 'L15_conn', 'L15_conn.feb')
PRESSURES = ('Load-AVW', 'Load-PVW', 'Load-PeB-top', 'Load-top', 'Load-LA')
NEAR_ZERO = ('PM_PeB_Left_', 'PM_PeB_Right_', 'PM_avw_bottom_left_', 'PM_avw_bottom_right_', 'PeB-constrin')
WANT = set(sys.argv[1:])


def want(name):
    return not WANT or name in WANT


def fresh(name, base=BASE):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(base)


def src(name):
    """As build_batch29.src: L11_stab with Load-LA back at 0.014 MPa, velocity plotted."""
    m = fresh(name)
    for sl in m.root.find('Loads'):
        if sl.get('name') in PRESSURES:
            p = sl.find('pressure')
            if abs(float(p.text) - 0.014) > 1e-9:
                old = p.text
                p.text = '0.014'
                m.log.append(f'{sl.get("name")} pressure {old} -> 0.014 MPa: back at the Abaqus value (all five '
                             f'pressures are now the source loads)')
    assert all(abs(float(s.find('pressure').text) - 0.014) < 1e-9 for s in m.root.find('Loads')
               if s.get('name') in PRESSURES)
    plot = m.root.find('Output').find('plotfile')
    if not any(v.get('type') == 'velocity' for v in plot):
        ET.SubElement(plot, 'var', type='velocity')
        m.log.append('output only: velocity added to the plot file (no model change)')
    return m


def cutback(m):
    m.set_solver(aggressiveness=1, cutback=0.5)
    m.log.append('solver only: a real cutback (FEBio reads cutback only with aggressiveness 1): each retry halves the '
                 'step instead of taking dt0/21 off it')


def parcus_to_connectors(m):
    doms = [d.get('name') for d in m.root.find('MeshDomains')]
    lofts = [d for d in ('P-arcus-L_fan', 'P-arcus-R_fan') if d in doms]
    assert len(lofts) == 2, doms
    m.remove_domains(lofts)
    m.log.append(f'removed the P-arcus lofts {lofts} (merged with the arcus chain nodes; no contacts on them)')
    m.add_parcus_connectors(extend='constant')


def unsplit_chain(m):
    """Undo parcus_edge_to_chain (batch 19): the 8 + 8 P-arcus loft edge nodes it put into the posterior-arcus spring
    chains sit between two collinear springs with no element and no mass; without the loft they are singular (L16_allconn
    failed its first step: 850 negative jacobians at ~zero load). Each such node's two springs -> one spring between its
    neighbours, as the source's T3D2 element (strain-measure law, so the same law)."""
    mesh = m.root.find('Mesh')
    elem_nodes = {int(v) for blk in mesh.findall('Elements') for e in blk for v in e.text.split(',')}
    for side in ('Left', 'Right'):
        ds = next(d for d in mesh.findall('DiscreteSet') if d.get('name') == f'Posterior_Arcus_{side}_springs')
        merged = []
        while True:
            pairs = [(e, tuple(int(v) for v in e.text.split(','))) for e in ds.findall('delem')]
            deg = {}
            for e, p in pairs:
                for n in p:
                    deg.setdefault(n, []).append((e, p))
            lone = sorted(n for n, es in deg.items() if len(es) == 2 and n not in elem_nodes)
            if not lone:
                break
            x = lone[0]
            (e1, p1), (e2, p2) = deg[x]
            a = p1[0] if p1[1] == x else p1[1]
            b = p2[1] if p2[0] == x else p2[0]
            e1.text = f'{a},{b}' if p1[1] == x else f'{b},{a}'
            ds.remove(e2)
            merged.append(x)
        n = len(ds.findall('delem'))
        m.log.append(f'Posterior_Arcus_{side}_springs: the {len(merged)} loft-edge nodes {merged} taken out of the chain '
                     f'again (each one\'s two springs -> one, the source element; strain-measure law unchanged): '
                     f'{n + len(merged)} -> {n} springs')


if want('L16_nzconn'):
    m = src('L16_nzconn'); m.lofts_to_connectors(NEAR_ZERO); emit('L16_nzconn', m)
if want('L16_nzconn_aggr'):
    m = src('L16_nzconn_aggr'); m.lofts_to_connectors(NEAR_ZERO); cutback(m); emit('L16_nzconn_aggr', m)
if want('L16_aggr'):
    m = src('L16_aggr'); cutback(m); emit('L16_aggr', m)
if want('L16_allconn'):
    m = fresh('L16_allconn', CONN); parcus_to_connectors(m); emit('L16_allconn', m)
if want('L16_allconn2'):
    m = fresh('L16_allconn2', CONN); parcus_to_connectors(m); unsplit_chain(m); emit('L16_allconn2', m)
if want('L16_conn_aggr'):
    m = fresh('L16_conn_aggr', CONN); cutback(m); emit('L16_conn_aggr', m)
