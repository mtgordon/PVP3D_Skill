"""What turns the perineal body (_PickedSet66) about x: force and moment of each thing attached to it, about its current
centroid, at a few converged states.

usage: py -3.10 peb_moments.py RUN [RUN ...] [--t 0.5,0.66,0.8,1.0] [--set _PickedSet66]
Groups: every discrete set with an end on the body (force from its table at the current elongation), the pressure loads
on the body's facets (p x load curve), and the contact pressure on the body's facets of every plotted contact surface.
M_-x = the moment about -x in N mm: positive drives the user's sense (cranial end forward and down), negative resists it.
The rest (the PVW through the 108 shared nodes, inertia) is the remainder: minus the sum, if the body is in equilibrium.
"""
import os
import struct
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RUNS_DIR  # noqa: E402
from febmodel import Feb  # noqa: E402
from xplt import Xplt, chunks  # noqa: E402

PLT_MESH = 0x01040000
PLT_SURFACE_SECTION = 0x01043000
PLT_SURFACE = 0x01043100
PLT_SURFACE_HDR = 0x01043101
PLT_SURFACE_ID = 0x01043102
PLT_SURFACE_FACES = 0x01043103
PLT_SURFACE_NAME = 0x01043104
PLT_SURFACE_MAX_FACET_NODES = 0x01043105
PLT_FACE_LIST = 0x01043200
PLT_FACE = 0x01043201


def read_surfaces(path):
    """[(id, name, faces (nf, 2 + maxn) int array: [face id, nnodes, 0-based node indices...])] from the mesh chunk."""
    out = []
    with open(path, 'rb') as f:
        f.read(4)
        cid, size = struct.unpack('<II', f.read(8))
        f.seek(size, 1)                       # root
        cid, size = struct.unpack('<II', f.read(8))
        assert cid == PLT_MESH, hex(cid)
        buf = f.read(size)
    for cid, o, s in chunks(buf):
        if cid != PLT_SURFACE_SECTION:
            continue
        for c2, o2, s2 in chunks(buf, o, o + s):
            if c2 != PLT_SURFACE:
                continue
            sid = name = None
            maxn = 4
            faces = []
            for c3, o3, s3 in chunks(buf, o2, o2 + s2):
                if c3 == PLT_SURFACE_HDR:
                    for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                        if c4 == PLT_SURFACE_ID:
                            sid = struct.unpack_from('<I', buf, o4)[0]
                        elif c4 == PLT_SURFACE_MAX_FACET_NODES:
                            maxn = struct.unpack_from('<I', buf, o4)[0]
                        elif c4 == PLT_SURFACE_NAME:
                            raw = buf[o4:o4 + s4]
                            n = struct.unpack_from('<I', raw, 0)[0] if s4 >= 4 else 0
                            name = (raw[4:4 + n] if 0 < n <= s4 - 4 else raw).split(b'\0')[0].decode(errors='replace')
                elif c3 == PLT_FACE_LIST:
                    for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                        if c4 == PLT_FACE:
                            faces.append(np.frombuffer(buf, dtype=np.int32, count=s4 // 4, offset=o4).copy())
            out.append((sid, name, np.array(faces)))
    return out


def area_vectors(xyz, faces):
    """Area vector (outward by node order) and centroid of each face row [id, nn, nodes...]."""
    av = np.zeros((len(faces), 3))
    cen = np.zeros((len(faces), 3))
    for i, fr in enumerate(faces):
        nn = fr[1]
        p = xyz[fr[2:2 + nn]]
        if nn == 3:
            av[i] = 0.5 * np.cross(p[1] - p[0], p[2] - p[0])
        else:
            av[i] = 0.5 * np.cross(p[2] - p[0], p[3] - p[1])
        cen[i] = p.mean(axis=0)
    return av, cen


def smooth_step(t, t0=0.0, t1=1.0):
    s = np.clip((t - t0) / (t1 - t0), 0, 1)
    return s * s * (3 - 2 * s)


def main():
    args = sys.argv[1:]
    times = (0.5, 0.66, 0.8, 1.0)
    dom = '_PickedSet66'
    if '--t' in args:
        i = args.index('--t')
        times = tuple(float(v) for v in args[i + 1].split(','))
        args = args[:i] + args[i + 2:]
    if '--set' in args:
        i = args.index('--set')
        dom = args[i + 1]
        args = args[:i] + args[i + 2:]
    for run in args:
        febp = os.path.join(RUNS_DIR, run, run + '.feb')
        xp = os.path.join(RUNS_DIR, run, run + '.xplt')
        fe = Feb(febp)
        root = fe.root
        x = Xplt(xp)
        idx = {int(n): i for i, n in enumerate(x.node_ids)}
        body_ids = set(fe.domain_nodes(dom))
        body = np.array(sorted(idx[n] for n in body_ids))
        bodyset = set(body.tolist())
        # discrete sets touching the body
        dmats = root.find('Discrete').findall('discrete_material') if root.find('Discrete') is not None else []
        tables = []
        for dm in dmats:
            pts = np.array([[float(v) for v in p.text.split(',')] for p in dm.iter('pt')])
            tables.append((dm.get('name'), float(dm.findtext('scale') or 1), (dm.findtext('measure') or 'elongation').strip(), pts))
        groups = []
        for d in (root.find('Discrete').findall('discrete') if root.find('Discrete') is not None else []):
            pairs = [(idx[a], idx[b]) for a, b in fe.discsets.get(d.get('discrete_set'), [])]
            on = [(a, b) if a in bodyset else (b, a) for a, b in pairs if (a in bodyset) != (b in bodyset)]
            if on:
                groups.append((d.get('discrete_set'), np.array(on), tables[int(d.get('dmat')) - 1]))
        # pressure loads on the body
        loads = []
        lc_scale = {}
        ld = root.find('LoadData')
        for lcn in ld.findall('load_controller') if ld is not None else []:
            pts = [tuple(float(v) for v in p.text.split(',')) for p in lcn.iter('pt')]
            lc_scale[lcn.get('id')] = pts
        for sl in root.find('Loads').findall('surface_load'):
            if sl.get('type') != 'pressure':
                continue
            facets = [c for _, c in fe.surfaces.get(sl.get('surface'), []) if all(n in body_ids for n in c)]
            if not facets:
                continue
            pe = sl.find('pressure')
            p0 = float(pe.text)
            pts = lc_scale.get(pe.get('lc'))
            fr = np.array([[0, len(c)] + [idx[n] for n in c] + [0] * (4 - len(c)) for c in facets])
            loads.append((sl.get('name'), p0, pts, fr))
        # contact surfaces in the plot file: facets fully on the body
        surfs = []
        for sid, name, faces in read_surfaces(xp):
            if len(faces) == 0:
                continue
            on = np.array([all(int(n) in bodyset for n in fr[2:2 + fr[1]]) for fr in faces])
            if on.any():
                surfs.append((sid, name, faces, on))
        conv = [k for k, st in enumerate(x.states) if st[1] == 0]
        tt = np.array([x.states[k][0] for k in conv])
        print(f'== {run}: {dom} ({len(body)} nodes); M_-x > 0 drives cranial-forward-down, < 0 resists  [N, N mm]')
        ks = sorted({conv[int(np.argmin(np.abs(tt - t)))] for t in times if t <= tt[-1] + 1e-6} | {conv[-1]})
        for k in ks:
            t = x.states[k][0]
            u = x.var(k, 'displacement')
            xyz = x.X + u
            c = xyz[body].mean(axis=0)
            rows = []
            for name, on, (mname, scale, measure, pts) in groups:
                a, b = on[:, 0], on[:, 1]
                L0 = np.linalg.norm(x.X[a] - x.X[b], axis=1)
                dv = xyz[b] - xyz[a]
                L = np.linalg.norm(dv, axis=1)
                e = L - L0 if measure == 'elongation' else (L - L0) / L0
                fmag = scale * np.interp(e, pts[:, 0], pts[:, 1])      # constant extension both ends
                F = (fmag / np.maximum(L, 1e-12))[:, None] * dv
                M = np.cross(xyz[a] - c, F)
                rows.append((f'{name} ({len(on)})', F.sum(axis=0), M.sum(axis=0),
                             f'elong med {np.median(e):+.2f} max {e.max():+.2f}; f max {fmag.max():.3f}'))
            for name, p0, pts, fr in loads:
                if pts:
                    tp = np.array(pts)
                    lf = smooth_step(t, tp[0, 0], tp[-1, 0]) * (tp[-1, 1] - tp[0, 1]) + tp[0, 1]
                else:
                    lf = 1.0
                av, cen = area_vectors(xyz, fr)
                F = -p0 * lf * av
                M = np.cross(cen - c, F)
                rows.append((f'{name} (p {p0 * lf:.4f}, {len(fr)} facets)', F.sum(axis=0), M.sum(axis=0), ''))
            cp = x.var(k, 'contact pressure')
            for sid, name, faces, on in surfs:
                if not isinstance(cp, dict) or sid not in cp:      # keyed by the 1-based surface id
                    continue
                arr = np.asarray(cp[sid]).ravel()
                if arr.size != len(faces):
                    rows.append((f'contact {name}', np.zeros(3), np.zeros(3), f'(data size {arr.size} != {len(faces)} faces)'))
                    continue
                av, cen = area_vectors(xyz, faces[on])
                F = -arr[on][:, None] * av
                M = np.cross(cen - c, F)
                rows.append((f'contact {name} ({int((arr[on] > 0).sum())}/{on.sum()} facets touching)', F.sum(axis=0),
                             M.sum(axis=0), f'p max {arr[on].max():.4f}'))
            tot_F = sum(r[1] for r in rows)
            tot_M = sum(r[2] for r in rows)
            print(f'-- t {t:.4f}: centroid ({c[0]:+.1f}, {c[1]:+.1f}, {c[2]:+.1f})')
            for name, F, M, note in rows:
                print(f'   {name:52s} F (y, z) ({F[1]:+7.2f}, {F[2]:+7.2f})  M_-x {-M[0]:+8.1f}   {note}')
            print(f'   {"sum of the above":52s} F (y, z) ({tot_F[1]:+7.2f}, {tot_F[2]:+7.2f})  M_-x {-tot_M[0]:+8.1f}   '
                  f'(the PVW and the rest carry the opposite)')


if __name__ == '__main__':
    main()
