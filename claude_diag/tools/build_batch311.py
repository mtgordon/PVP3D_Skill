"""Batch 311 (2026-10-17; the user saw the curves wrinkle into folds where their distal end bends ~90 deg near the PeB (the bend
itself "looks okay") and said "try both": (1) a pre-bend of the curve end, (2) a fold contact (batch 310). This is (1)).
NOT IN SOURCE (rest geometry change): the distal canal end - every node of the walls (AVW, PVW), the PeB and the curves
(canal_wrap_hex_*) with rest z below Z0 - is moved by F x its SHAPE change in REF_RUN at t TREF: the deformed positions of the
nodes with rest z < -15 are mapped back by their best-fit rigid motion (Kabsch: the descent and tilt removed), and the rest of
the displacement (the bend) is scaled by F and blended in with a smoothstep from rest z Z0 (0) to Z1 (1). Nodes in a fixed BC
set do not move. So the curve end starts part-way bent and needs less bending (and less compression on the inside of the bend)
to reach the same shape. Springs keep their own rest lengths; ties follow because the field comes from a solution that
satisfied them.
--all: also move the other nodes (lofts, chains) within 10 mm of a moved canal node by that node's move, fading to 0 at 10 mm
(PM arch, levator and fixed nodes never; L311 moved only walls / PeB / curves and flipped fan-loft triangles).
usage: py -3.10 build_batch311.py SRC_RUN NEW_RUN REF_RUN [--t 0.6] [--f 0.5] [--z0 -12] [--z1 -22] [--all]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402
from xplt_reader import Xplt, converged_state_indices  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
DOMS = ('_PickedSet347', '_PickedSet64', '_PickedSet66', 'canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')


def main(src, new, ref, tref, f, z0, z1, allnodes=False):
    t, X, E = load(src)
    canal = sorted({v for nm, c in E.values() if nm in DOMS for v in c})
    nodes = canal
    if allnodes:   # every node in the region except the PM arch, the levator (LA*) and fixed nodes: lofts / chains follow
        skip = {v for nm, c in E.values() if nm == 'PM_Plane' or nm.startswith('LA') for v in c}
        nodes = sorted(n for n in X if n not in skip)
    fixed = set()
    for name in re.findall(r'<bc name="[^"]*" node_set="([^"]+)" type="zero displacement">', t):
        m = re.search(rf'<NodeSet name="{re.escape(name)}">([^<]*)</NodeSet>', t)
        if m:
            fixed |= {int(v) for v in re.findall(r'\d+', m.group(1))}
    x = Xplt(os.path.join(RUNS, ref, ref + '.xplt'))
    idx = converged_state_indices(x, os.path.join(RUNS, ref, ref + '.log'))
    tt = np.array([x.states[i][0] for i in idx]); s = idx[int(np.argmin(abs(tt - tref)))]
    pos = {int(n): k for k, n in enumerate(x.node_ids)}
    U = x.var(s, 'displacement')
    reg = [n for n in canal if X[n][2] < -15 and n in pos]   # the rigid fit always from the canal end (walls / PeB / curves)
    P = np.array([X[n] for n in reg]); Q = P + np.array([U[pos[n]] for n in reg])
    cp, cq = P.mean(0), Q.mean(0)
    Hm = (P - cp).T @ (Q - cq); u_, _, vt = np.linalg.svd(Hm); d = np.sign(np.linalg.det(vt.T @ u_.T))
    Rm = vt.T @ np.diag([1, 1, d]) @ u_.T            # Q ~ Rm P + (cq - Rm cp)
    moved, dmax = 0, 0.0
    newX = {}
    for n in canal:
        z = X[n][2]
        if z >= z0 or n in fixed or n not in pos:
            continue
        w = min(1.0, (z0 - z) / (z0 - z1)); w = w * w * (3 - 2 * w)
        q = X[n] + U[pos[n]]
        back = Rm.T @ (q - (cq - Rm @ cp))           # deformed shape mapped back by the rigid fit
        dv = f * w * (back - X[n])
        newX[n] = X[n] + dv; moved += 1; dmax = max(dmax, np.linalg.norm(dv))
    if allnodes:   # other nodes (lofts, chains, anchors of patches...): the move of the nearest canal node, fading to 0 over DFADE mm
        from scipy.spatial import cKDTree
        cn = list(newX); tree = cKDTree(np.array([X[n] for n in cn])); cs = set(canal)
        for n in nodes:
            if n in cs or n in fixed:
                continue
            dist, k = tree.query(X[n])
            if dist >= 10.0:
                continue
            w = 1 - dist / 10.0; w = w * w * (3 - 2 * w)
            dv = w * (newX[cn[k]] - X[cn[k]])
            if np.linalg.norm(dv) > 1e-9:
                newX[n] = X[n] + dv; moved += 1
    def rep(m):
        n = int(m.group(1))
        if n in newX:
            p = newX[n]; return f'<node id="{n}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>'
        return m.group(0)
    t = re.sub(r'<node id="(\d+)">[^<]*</node>', rep, t)
    ang = np.degrees(np.arccos(np.clip((np.trace(Rm) - 1) / 2, -1, 1)))
    msg = (f'NOT IN SOURCE (rest geometry): pre-bend of the distal canal end from {ref} at t {x.states[s][0]:.3f}: rigid motion removed '
           f'(best-fit rotation {ang:.1f} deg + translation {np.linalg.norm(cq - cp):.1f} mm), the shape change x{f:g}, smoothstep from '
           f'rest z {z0:g} to {z1:g}; {moved} nodes moved (' + ('canal + nearby lofts / chains, fading over 10 mm' if allnodes else 'walls / PeB / curves') + f') (max {dmax:.2f} mm); fixed-BC nodes kept')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch311.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], a[2], g('--t', 0.6), g('--f', 0.5), g('--z0', -12.0), g('--z1', -22.0), '--all' in a)
