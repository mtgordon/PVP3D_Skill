"""Shell pinch screen using FEBio 4's default shell geometry (FEBio User Manual 3.6.2.2):
the nodes define the shell TOP face and the bottom face lies a full thickness t away along the
negative (averaged) nodal normal. Through-thickness Gauss points (2-point rule) sit at
t*(1 -+ 1/sqrt(3))/2 below the nodes.

The skill's shell_pinch_check.py offsets +-t/2 about the nodes (the pre-2.6 mid-surface
convention, which FEBio 4 uses only for elastic-shell-old), so it understates the pinch by half.

For each element: signed area of the offset surface / area of the node surface, at the
bottom face and at the two Gauss-point levels. 1 = prism-like, < 0.25 badly pinched, <= 0 inverted.
Domains whose ShellDomain already has shell_normal_nodal=0 use element normals (never pinched).

usage: py -3.10 shell_pinch_top.py model.feb [--domains A,B] [--convention top|mid] [--list N]
"""
import argparse
from collections import defaultdict

import numpy as np

from febmodel import Feb


def elem_normal(P):
    n = np.cross(P[1] - P[0], P[2] - P[0]) if len(P) == 3 else np.cross(P[2] - P[0], P[3] - P[1])
    return n / np.linalg.norm(n)


def signed_area(P, nref):
    if len(P) == 3:
        return 0.5 * np.cross(P[1] - P[0], P[2] - P[0]) @ nref
    return 0.5 * np.cross(P[2] - P[0], P[3] - P[1]) @ nref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('feb')
    ap.add_argument('--domains')
    ap.add_argument('--convention', choices=('top', 'mid'), default='top')
    ap.add_argument('--list', type=int, default=10)
    a = ap.parse_args()
    f = Feb(a.feb)
    thick, nodal = {}, {}
    for dom in f.section('MeshDomains'):
        if dom.tag == 'ShellDomain':
            thick[dom.get('name')] = float(dom.find('shell_thickness').text.split(',')[0])
            sn = dom.find('shell_normal_nodal')
            nodal[dom.get('name')] = sn is None or sn.text.strip() != '0'
    names = a.domains.split(',') if a.domains else sorted(thick)
    # nodal normals: average over ALL shell elements sharing the node (manual wording)
    acc = defaultdict(lambda: np.zeros(3))
    for n in thick:
        for conn in f.elem_blocks[n][1].values():
            nn = elem_normal(np.array([f.nodes[i] for i in conn]))
            for i in conn:
                acc[i] += nn
    if a.convention == 'top':
        levels = {'bottom': 1.0, 'gp_low': (1 + 1 / np.sqrt(3)) / 2, 'gp_high': (1 - 1 / np.sqrt(3)) / 2}
    else:
        levels = {'bottom': 0.5, 'gp_low': 0.5 / np.sqrt(3), 'top': -0.5, 'gp_high': -0.5 / np.sqrt(3)}
    print(f'convention: {a.convention} (offset below the nodes, in units of t: {levels})')
    worst = []
    for n in names:
        t = thick[n]
        rows = []
        for eid, conn in f.elem_blocks[n][1].items():
            P = np.array([f.nodes[i] for i in conn])
            ne = elem_normal(P)
            D = np.array([acc[i] for i in conn]) if nodal[n] else np.tile(ne, (len(conn), 1))
            D /= np.linalg.norm(D, axis=1)[:, None]
            a0 = signed_area(P, ne)
            r = {k: signed_area(P - s * t * D, ne) / a0 for k, s in levels.items()}
            face = min(v for k, v in r.items() if not k.startswith('gp'))
            gp = min(v for k, v in r.items() if k.startswith('gp'))
            rows.append((face, gp, eid))
        R = np.array([(x[0], x[1]) for x in rows])
        print(f'{n:26s} t={t:<5g} nodal_normals={int(nodal[n])}  face<=0: {int((R[:, 0] <= 0).sum()):4d}  '
              f'face<0.25: {int((R[:, 0] < 0.25).sum()):4d}  GP<=0: {int((R[:, 1] <= 0).sum()):4d}  '
              f'min face {R[:, 0].min():7.3f}  min GP {R[:, 1].min():7.3f}')
        worst += [(x[0], x[1], n, x[2]) for x in rows]
    print('worst elements (face ratio, GP ratio, domain, id):')
    for w in sorted(worst)[:a.list]:
        print(f'   {w[0]:7.3f} {w[1]:7.3f}  {w[2]:24s} {w[3]}')


if __name__ == '__main__':
    main()
