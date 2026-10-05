"""Perineal-body (_PickedSet66) rigid rotation: the best-fit (Kabsch) rotation of its nodes from the reference to each
converged state, as a rotation vector, plus the centroid translation.

usage: py -3.10 peb_rotation.py RUN [RUN ...] [--t 0.25,0.5,0.75,1.0] [--set NAME[,NAME]] [--all]
  rot_-x = the rotation about -x in degrees: the user's sense (2026-09-26), cranial end forward (+y) and down (-z);
           a node above the centroid moves +y, a node in front of it moves -z. Positive = the wanted sense.
  rotvec = the full rotation vector (deg about x, y, z; right-hand rule).
  fit    = rms distance (mm) of the deformed nodes from the rigidly moved reference: how far from rigid the body moved.
  --all  every converged state (else the states nearest the --t times, plus the last converged one).
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from xplt import Xplt  # noqa: E402


def kabsch(P, Q):
    """Rotation R minimising |R (P - Pc) - (Q - Qc)|; returns R, rotation vector (rad), rms residual."""
    Pc, Qc = P.mean(axis=0), Q.mean(axis=0)
    A, B = P - Pc, Q - Qc
    H = A.T @ B
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    ang = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    if ang < 1e-9:
        rv = np.zeros(3)
    else:
        axis = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / (2 * np.sin(ang))
        rv = axis * ang
    res = np.sqrt(((A @ R.T - B) ** 2).sum(axis=1).mean())
    return R, rv, res


def region_nodes(run, x, sets):
    names = [d.get('name') for d in ET.parse(os.path.join(RUNS_DIR, run, run + '.feb')).getroot().find('MeshDomains')]
    s = set()
    for i, d in enumerate(x.domains):
        if i < len(names) and names[i] in sets:
            s.update(int(n) for c in d['conn'] for n in c)
    return np.array(sorted(s))


def main():
    args = sys.argv[1:]
    times = (0.25, 0.5, 0.75, 1.0)
    sets = ('_PickedSet66',)
    every = False
    if '--t' in args:
        i = args.index('--t')
        times = tuple(float(v) for v in args[i + 1].split(','))
        args = args[:i] + args[i + 2:]
    if '--set' in args:
        i = args.index('--set')
        sets = tuple(args[i + 1].split(','))
        args = args[:i] + args[i + 2:]
    if '--all' in args:
        every = True
        args.remove('--all')
    for run in args:
        x = Xplt(os.path.join(RUNS_DIR, run, run + '.xplt'))
        nodes = region_nodes(run, x, sets)
        conv = [k for k, st in enumerate(x.states) if st[1] == 0]
        if not conv:
            print(f'== {run}: no converged state')
            continue
        tt = np.array([x.states[k][0] for k in conv])
        P = x.X[nodes]
        print(f'== {run} ({"+".join(sets)}, {len(nodes)} nodes; last converged t = {tt[-1]:.4f})')
        print(f'{"t":>7} {"rot_-x":>7} | rotvec x y z [deg]      | centroid u x y z [mm]   | fit rms [mm]')
        if every:
            ks = conv
        else:
            ks = []
            for t in times:
                if t <= tt[-1] + 1e-6:
                    ks.append(conv[int(np.argmin(np.abs(tt - t)))])
            ks.append(conv[-1])
            ks = sorted(set(ks))
        for k in ks:
            u = x.var(k, 'displacement')[nodes]
            R, rv, res = kabsch(P, P + u)
            rv = np.degrees(rv)
            cu = u.mean(axis=0)
            print(f'{x.states[k][0]:7.4f} {-rv[0]:+7.2f} | ({rv[0]:+6.2f}, {rv[1]:+6.2f}, {rv[2]:+6.2f}) '
                  f'| ({cu[0]:+5.1f}, {cu[1]:+5.1f}, {cu[2]:+5.1f}) | {res:5.2f}')


if __name__ == '__main__':
    main()
