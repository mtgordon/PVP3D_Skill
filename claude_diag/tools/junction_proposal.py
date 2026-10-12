"""The proposed conforming tube junction (2026-10-10, for the user's picture; NOT a model yet). On each side the curve's last two
PVW-side hex rings (12 hexes each; they run through the slanted lateral face of the PVW's last layer and are held only by the
canal_wrap_embed point ties) are re-placed so that:
  - the curve's PVW-face column (the column whose nodes the PVW shares higher up) runs down ALONG the PVW's slanted lateral face
    and ends ON the PVW / PeB corner nodes of the PVW's distal face (shared node IDs: PVW and PeB already share them);
  - the middle node ring lands half-way down that face (new nodes there: the PVW's last layer split in two through its depth
    so that ring has PVW nodes to share, i.e. the PVW end remeshed in that one layer);
  - the other columns of the U follow that shift, tapering to nothing at the AVW / PVW midpoint (column 3), so the U keeps its
    shape and no curve hex sits inside the PVW / PeB any more; canal_wrap_embed removed.
Hanging curve nodes on the PVW face (the curve has 5 rows through the thickness, the PVW face 4) stay on PVW edges as now
(canal_wrap_ends constraints). Writes the proposal's JSON for the viewer and reports the overlap check.
usage: py -3.10 junction_proposal.py RUN OUT.json"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load, frame, overlap  # noqa: E402
from apex_close import HF  # noqa: E402

RIGHT = [(27795, 27806), (27807, 27818), (27819, 27830)]   # rings A (last one sharing PVW nodes), B, C (end), right side


def ring_faces(E, L, rng):
    els = [E[k][1] for k in range(rng[0], rng[1] + 1)]
    top, bot = set(), set()
    for c in els:
        a, b = np.mean([L[v][2] for v in c[:4]]), np.mean([L[v][2] for v in c[4:]])
        (top if a > b else bot).update(c[:4]); (bot if a > b else top).update(c[4:])
    return top, bot


def mirror_ids(E, L, ids):
    """the left-side element of each right-side element (mirrored centroid)"""
    cand = [(k, np.mean([L[v] for v in c], 0)) for k, (n, c) in E.items() if n == 'canal_wrap_hex_pvw']
    K = np.array([k for k, _ in cand]); C = np.array([p for _, p in cand])
    out = []
    for k in ids:
        p = np.mean([L[v] for v in E[k][1]], 0) * [-1, 1, 1]
        out.append(int(K[np.argmin(np.linalg.norm(C - p, axis=1))]))
    return out


def grid(top, mid, end, L, pvw):
    """columns of the ring by distance from the PVW face: col 0 = the PVW-face column (the nodes of `top` shared with the PVW)"""
    shared = [v for v in top if v in pvw]
    face_x = np.mean([abs(L[v][0]) for v in shared])
    def cols(S):
        xs = sorted(S, key=lambda v: abs(L[v][0]) - face_x + 0.0 * L[v][1])
        return xs
    return shared, face_x


def main(run, out):
    t, X, E = load(run)
    c0, B = frame(X, E)
    L = {k: (p - c0) @ B for k, p in X.items()}
    dn = lambda d: {v for n, c in E.values() if n == d for v in c}   # noqa: E731
    pvw, peb = dn('_PickedSet64'), dn('_PickedSet66')
    newL = dict(L)
    report = []
    shared_all, new_nodes = [], []
    for side, rings in ((1, RIGHT), (-1, None)):
        if rings is None:
            rings = []
            for a, b in RIGHT:
                m = mirror_ids(E, L, range(a, b + 1)); rings.append((min(m), max(m)))
        topA, botA = ring_faces(E, L, rings[0])
        top, mid0 = ring_faces(E, L, rings[1])
        _, end0 = ring_faces(E, L, rings[2])
        mid0, end0 = sorted(mid0), sorted(end0)
        face_top = sorted(top & pvw, key=lambda v: L[v][1])              # the 4 shared PVW-face nodes at depth ~5.3
        # the PVW / PeB corner under each shared face node: the most lateral PVW node on the distal face in its AP row
        corner = []
        for v in face_top:
            cand = [u for u in pvw if L[u][0] * side > 3 and L[u][2] < L[v][2] - 1.5 and abs(L[u][1] - L[v][1]) < 0.45]
            corner.append(max(cand, key=lambda u: abs(L[u][0])))
        # curve end column = the end-ring nodes nearest the PVW face (5 rows), matched to the corner nodes by AP
        end_face = sorted(sorted(end0, key=lambda v: abs(L[v][0]))[:5], key=lambda v: L[v][1])
        mid_face = sorted(sorted(mid0, key=lambda v: abs(L[v][0]))[:5], key=lambda v: L[v][1])
        cAP = np.array([L[v][1] for v in corner])
        def on_face(v, frac):
            """point on the PVW's slanted lateral face at this node's AP, frac 0 = top (depth 5.3), 1 = distal corner"""
            ap = L[v][1]
            j = int(np.clip(np.searchsorted(cAP, ap) - 1, 0, len(cAP) - 2))
            w = np.clip((ap - cAP[j]) / (cAP[j + 1] - cAP[j]), 0, 1)
            topp = (1 - w) * L[face_top[j]] + w * L[face_top[j + 1]]
            botp = (1 - w) * L[corner[j]] + w * L[corner[j + 1]]
            return (1 - frac) * topp + frac * botp
        shift_end = {v: on_face(v, 1.0) - L[v] for v in end_face}
        shift_mid = {v: on_face(v, 0.5) - L[v] for v in mid_face}
        # taper the shift across the U: each node takes the shift of its nearest face-column node, scaled down linearly with
        # its distance from the face column, to nothing at the far end of the ring (the AVW / PVW midpoint, column 3)
        for ring_nodes, face, sh in ((end0, end_face, shift_end), (mid0, mid_face, shift_mid)):
            F = np.array([L[v] for v in face])
            near = {v: int(np.argmin(np.linalg.norm(F - L[v], axis=1))) for v in ring_nodes}
            d = {v: np.linalg.norm(F[near[v]] - L[v]) for v in ring_nodes}
            D = max(d.values())
            for v in ring_nodes:
                newL[v] = L[v] + (1 - d[v] / D) * sh[face[near[v]]]
        shared_all += [v for v in end_face]
        new_nodes += [v for v in mid_face]
        report.append(f'side {side:+d}: rings {rings}; PVW face top {face_top}; PVW / PeB corner nodes {corner} '
                      f'(PeB shares {sum(c in peb for c in corner)} of them); end-face shift '
                      f'{np.round(np.mean([np.linalg.norm(s) for s in shift_end.values()]), 2)} mm mean')
    # overlap check on the proposed geometry
    Xn = {k: c0 + B @ p for k, p in newL.items()}
    before, after = overlap(X, E), overlap(Xn, E)
    report.append(f'curve hexes with the centroid inside a PVW / PeB hex: now {len(before)} -> proposal {len(after)} {sorted(after)}')
    print('\n'.join(report))
    json.dump({'newL': {int(k): [round(float(c), 3) for c in v] for k, v in newL.items() if np.linalg.norm(v - L[k]) > 1e-9},
               'shared': shared_all, 'mid': new_nodes, 'report': report}, open(out, 'w'))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
