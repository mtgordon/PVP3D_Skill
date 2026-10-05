"""Batch 143 (2026-10-01, the user: "It looks like there's a flare at the end that you added. Can you make it just continue the
curved edges as before and have it intersect with the pvw and peb?").
As L142_thin275_narrow3_wrapend (the curves carried to the last seam station at the perineal-body end; 2.75 mm walls narrowed
3 mm; the curves' inner faces in the lumen contact; walls_LA + PVW_LA; seg_up 2; the PM as a structure; the 27 PM_conn-
family springs on the PM's inner arc), but:
  * no flare on the AVW, the cervix or the curves: they stay narrowed to the end, so each curve keeps its shape;
  * the PVW alone flares back to full width over the last 3.4 mm before the perineal-body junction (as L137), so the
    perineal body is not narrowed;
  * the curves are built on the walls with the PVW unflared (geometry A), then put in the final geometry (B: PVW flared):
    wherever a curve node or a curve end that followed a PVW column would be dragged by the PVW's flare, it keeps its
    geometry-A position and is embedded in the PVW / perineal-body hex it lies in (or the nearest one, extrapolated):
    u - sum N_i(xi) u_i = 0 per dof (constraint canal_wrap_embed, penalty 10, maxaug 10, as canal_wrap_ends). So the curves
    intersect the PVW and the perineal body there, joined to them, instead of flaring. Curve faces with an embedded node
    are left out of the lumen contact (they start inside the PVW).
  L143_thin275_narrow3_wrapend_noflare   built, then run
usage: py -3.10 build_batch143.py"""
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from build_batch90 import seam_nodes  # noqa: E402
from build_batch118 import fold  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM  # noqa: E402
from build_batch132 import seg_up2  # noqa: E402
from febmodel import Feb  # noqa: E402
import tube  # noqa: E402
import wall_reshape  # noqa: E402
import build_batch135 as b135  # noqa: E402
import build_batch136 as b136  # noqa: E402
import build_batch141 as b141  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
NAME = 'L143_thin275_narrow3_wrapend_noflare'
PRE, MID, AVG = NAME + '_tmppre', NAME + '_tmpmid', NAME + '_tmpavg'
TOOLS = os.path.dirname(os.path.abspath(__file__))
XI = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def shape(xi):
    return np.prod(1 + XI * xi, axis=1) / 8


def dshape(xi):
    d = np.zeros((8, 3))
    for k in range(3):
        o = [j for j in range(3) if j != k]
        d[:, k] = XI[:, k] * (1 + XI[:, o[0]] * xi[o[0]]) * (1 + XI[:, o[1]] * xi[o[1]]) / 8
    return d


def natural(P, p):
    xi = np.zeros(3)
    for _ in range(25):
        r = shape(xi) @ P - p
        J = P.T @ dshape(xi)
        try:
            dx = np.linalg.solve(J, r)
        except np.linalg.LinAlgError:
            break
        xi -= dx
        if np.linalg.norm(dx) < 1e-10:
            break
    return xi


def embed(p, hosts, HC):
    """(element node ids, weights, max |xi|) of the host hex containing p, or the nearest one (extrapolated)"""
    best = None
    for k in np.argsort(np.linalg.norm(HC - p, axis=1))[:12]:
        conn, P = hosts[k]
        xi = natural(P, p)
        m = np.abs(xi).max()
        if best is None or m < best[2]:
            best = (conn, shape(xi), m)
        if m <= 1.0 + 1e-6:
            break
    return best


if __name__ == '__main__':
    os.makedirs(SCRATCH, exist_ok=True)
    for nm in (NAME, PRE, MID, AVG):
        assert not os.path.exists(os.path.join(RUNS, nm)), f'{nm} exists; not overwriting'
    orig = os.path.join(RUNS, BASE, BASE + '.feb')
    wall_reshape.NARROW = 3.0
    spec = b136.flare_spec(orig)
    # geometry A: the walls narrowed, nothing flared (only to build the curves on)
    mA = Model7(orig)
    wall_reshape.reshape(mA)
    pA = os.path.join(SCRATCH, NAME + '_geomA.feb')
    mA.write(pA)
    # geometry B: as L137: the PVW flared back toward the perineal body, the PeB not narrowed
    mB = Model7(orig)
    repB = wall_reshape.reshape(mB, flare=spec, flare_doms=(tube.PVW,))
    XB = mB.nodes()
    m = Model113(pA)
    m.log.append(f'{b135.LABEL}: wall_reshape (narrowed 3 mm; the PVW alone flared back to full width over '
                 f'{[round(v[1], 1) for v in spec.values()]} mm before the PeB junction; the PeB not narrowed): {repB}')
    b135.NAME = NAME
    b135.EXTEND_DISTAL, b135.EXTENDED = True, {}
    jh = b135.wrap2(m, orig)
    XA = m.nodes()
    print(m.log[-1][:300])
    # put the walls / PeB in geometry B; the curves keep geometry A
    wall_reshape._set_coords(m, {n: XB[n] for n in XB})
    moved = {n for n in XB if np.linalg.norm(XB[n] - XA[n]) > 0.01}
    fB = Feb(orig)
    hosts = []
    for dm in (tube.PVW, tube.PEB):
        for c in fB.elem_blocks[dm][1].values():
            hosts.append((c, np.array([XB[n] for n in c])))
    HC = np.array([P.mean(0) for _, P in hosts])
    wrap_blocks = [b for b in m.mesh.findall('Elements') if b.get('name', '').startswith('canal_wrap_hex')]
    used = {int(v) for b in wrap_blocks for e in b for v in e.text.split(',')}
    dup_src = sorted(used & moved)                      # wall nodes the curves shared that the PVW flare moved
    nid0 = m.max_node_id() + 1
    dup_ids = m.add_nodes('canal_wrap_embed_nodes', [XA[n] for n in dup_src])
    dmap = dict(zip(dup_src, dup_ids))
    for b in wrap_blocks:
        for e in b:
            ids = [int(v) for v in e.text.split(',')]
            if any(i in dmap for i in ids):
                e.text = ','.join(str(dmap.get(i, i)) for i in ids)
    # canal_wrap_ends entries that follow a moved column: replaced by embedding
    cons = m.root.find('Constraints')
    ce = next(c for c in cons if c.get('name') == 'canal_wrap_ends')
    to_embed = set(dup_ids)
    drop = []
    for lc in ce.findall('linear_constraint'):
        nodes = [int(nd.get('id')) for nd in lc.findall('node')]
        if any(n in moved for n in nodes[1:]):
            to_embed.add(nodes[0])
            drop.append(lc)
    for lc in drop:
        ce.remove(lc)
    XA2 = {**XA, **{d: XA[s] for s, d in dmap.items()}}
    c = ET.SubElement(cons, 'constraint', name='canal_wrap_embed', type='linear constraint')
    ET.SubElement(c, 'tol').text = '0.01'
    ET.SubElement(c, 'penalty').text = '10'
    ET.SubElement(c, 'maxaug').text = '10'
    outside = []
    for n in sorted(to_embed):
        conn, w, mx = embed(XA2[n], hosts, HC)
        if mx > 1.0 + 1e-3:
            outside.append(round(float(mx), 2))
        for dof in ('x', 'y', 'z'):
            lc = ET.SubElement(c, 'linear_constraint')
            ET.SubElement(lc, 'node', id=str(n), bc=dof).text = '1'
            for q, wq in zip(conn, w):
                if abs(wq) > 1e-9:
                    ET.SubElement(lc, 'node', id=str(q), bc=dof).text = f'{-wq:.6f}'
    m.log.append(f'NOT IN SOURCE: the curves keep their unflared shape to the end and intersect the PVW / PeB: {len(dup_src)} '
                 f'wall nodes they shared were moved by the PVW flare and are replaced in the curves by duplicates at the '
                 f'unflared positions; {len(drop) // 3} curve end nodes that followed a moved column; all {len(to_embed)} '
                 f'embedded in the PVW / PeB hex they lie in (constraint canal_wrap_embed, penalty 10, maxaug 10); '
                 f'{len(outside)} lie outside every host hex (nearest one, extrapolated; max |xi| {max(outside) if outside else 0})')
    print(m.log[-1])
    seg_up2(m)
    step2(INNER + OUTER_BOTTOM)(m)
    inside = b135.walls_la(m)
    print(m.log[-1][:300])
    emit(PRE, m)
    b136.run_main(os.path.join(TOOLS, 'build_batch129.py'), {'BASE': PRE, 'NAME': MID})
    b136.run_main(os.path.join(TOOLS, 'build_batch129b.py'), {'SRC': MID, 'OLD': PRE, 'NAME': AVG})
    seam, _, _ = seam_nodes(Model7(orig))
    rp, _ = fold(Model7(orig))
    inv = {a: p for p, a in rp.items()}
    inv.update(b135.EXTENDED)
    src = os.path.join(RUNS, AVG, AVG + '.feb')
    f, t = Feb(src), open(src, encoding='utf-8').read()
    C = np.array([0.5 * (np.array(XA[a]) + np.array(XA[inv[a]])) for a in seam if a in inv])
    fp = [fc for fc in b141.inner_faces(f, 'canal_wrap_hex_pvw', C) if not any(n in to_embed for n in fc)]
    fa = [fc for fc in b141.inner_faces(f, 'canal_wrap_hex_avw', C) if not any(n in to_embed for n in fc)]
    t2 = b141.add_faces(b141.add_faces(t, 'SlidingElastic1Primary', fp), 'SlidingElastic1Secondary', fa)
    out = os.path.join(RUNS, NAME)
    os.makedirs(out)
    open(os.path.join(out, NAME + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
    notes = [open(os.path.join(RUNS, nm, nm + '.feb.changes.txt'), encoding='utf-8').read() for nm in (PRE, MID, AVG)]
    with open(os.path.join(out, NAME + '.feb.changes.txt'), 'w', encoding='utf-8') as fo:
        fo.write(f'base: {orig} (tools/build_batch143.py)\n' + '\n'.join(notes)
                 + f'\nthe curves\' inner faces in the lumen contact (faces with an embedded node left out): {len(fp)} '
                   f'PVW-half faces -> SlidingElastic1Primary, {len(fa)} AVW-half faces -> SlidingElastic1Secondary\n')
    for nm in (PRE, MID, AVG):
        shutil.move(os.path.join(RUNS, nm), os.path.join(SCRATCH, nm))
    print('inner faces', len(fp), len(fa), '; wall nodes inside the LA at rest', len(inside))
    print('wrote', os.path.join(out, NAME + '.feb'))
