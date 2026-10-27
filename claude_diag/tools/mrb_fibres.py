"""Muscle-fibre directions from a 3D Slicer scene (.mrb) onto the FEBio levator (2026-10-24, Stephen's turn; Mark's e-mail:
Michigan's "Perin Complex EMRLD 32318 New LA Fibers.mrb" holds levator fibre directions + the pelvic-floor segmentations).

An .mrb is a zip: a .mrml scene (XML listing every node: volumes, segmentations, models, markups, fibre bundles) plus the data
files (.nrrd / .seg.nrrd volumes, .vtk / .vtp / .ply / .stl models and fibre bundles, .mrk.json markups). Pure Python (numpy,
scipy only), so it runs anywhere the other tools run.

  inspect:  py -3.10 mrb_fibres.py inspect FILE.mrb
            lists the scene's nodes and every data file with what it holds (points / lines / vectors / volume size).
  map:      py -3.10 mrb_fibres.py map FILE.mrb RUN --fibres PATTERN --la-surface PATTERN [--out OUT.xml] [--max-dist 5]
            1. fibres: every polyline in the data files matching PATTERN (fibre bundles, curves) -> unit segment directions
               at the segment midpoints (Slicer RAS mm; LPS files are flipped to RAS);
            2. registration: the levator surface from the scene (model / segmentation-derived points matching
               --la-surface) is aligned to RUN's levator nodes (centroid + principal axes, then rigid ICP; the best of the
               8 axis-sign flips by ICP residual). Reported: RMS residual, scale ratio. NOTHING is assumed about the
               scanner frame vs the Abaqus frame, so check the residual (a good fit is ~1-3 mm for a 4 mm-thick shell);
            3. per levator shell element: the fibre segments within --max-dist of its centroid (after registration),
               weighted 1/d, their directions as an orientation tensor (sign-free: a fibre has no head or tail); the
               principal eigenvector, projected into the shell's tangent plane = the element's fibre direction; also the
               coherence (eigenvalue gap, 0..1) and the number of segments, so poorly covered elements are visible;
            4. writes <ElementData type="mat_axis" elem_set="LA_..."> blocks (a = fibre, d = the in-plane normal to it)
               for each levator domain, ready to paste into the .feb's <MeshData>, plus a CSV (element, centroid,
               direction, coherence, count). Elements with no fibre nearby get the nearest covered element's direction.
            The levator material must then be anisotropic for the directions to matter (now: isotropic Yeoh), e.g. a
            solid mixture of the Yeoh matrix + a tension-only fibre (fiber-exp-pow / fiber-NH) along mat_axis a: a
            separate, fitted change (NOT IN SOURCE: the source's levator is isotropic Marlow).
"""
import io
import json
import os
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile
import gzip
import zlib
import base64

import numpy as np
from scipy.spatial import cKDTree

LA_DOMAINS = ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri')


# ----------------------------------------------------------------------------------------------------------------- readers
def read_nrrd_header(b):
    end = b.find(b'\n\n')
    hdr = b[:end].decode('latin-1').splitlines()
    h = {}
    for line in hdr[1:]:
        if line.startswith('#') or ':' not in line:
            continue
        k, v = line.split(':', 1)
        h[k.strip().lower()] = v.lstrip('=').strip()
    return h, end + 2


def read_nrrd(b):
    h, off = read_nrrd_header(b)
    sizes = [int(v) for v in h['sizes'].split()]
    dt = {'unsigned char': 'u1', 'uchar': 'u1', 'uint8': 'u1', 'short': 'i2', 'int16': 'i2', 'unsigned short': 'u2',
          'ushort': 'u2', 'int': 'i4', 'int32': 'i4', 'unsigned int': 'u4', 'float': 'f4', 'double': 'f8',
          'signed char': 'i1', 'int8': 'i1'}[h['type']]
    endian = '<' if h.get('endian', 'little') == 'little' else '>'
    raw = b[off:]
    enc = h.get('encoding', 'raw')
    if enc in ('gzip', 'gz'):
        raw = gzip.decompress(raw)
    elif enc != 'raw':
        raise ValueError('nrrd encoding ' + enc)
    a = np.frombuffer(raw, dtype=endian + dt, count=int(np.prod(sizes))).reshape(sizes[::-1])
    return h, a


def _vtk_legacy(b):
    """VTK legacy POLYDATA (ascii or binary): POINTS and LINES (+ POLYGONS for surfaces)."""
    txt = b
    binary = b'\nBINARY' in b[:400]
    pos = 0
    out = {}

    def nextline():
        nonlocal pos
        e = txt.index(b'\n', pos)
        line = txt[pos:e].decode('latin-1').strip()
        pos = e + 1
        return line
    while pos < len(txt):
        try:
            line = nextline()
        except ValueError:
            break
        tok = line.split()
        if not tok:
            continue
        key = tok[0].upper()
        if key == 'POINTS':
            n, typ = int(tok[1]), tok[2].lower()
            if binary:
                dt = '>f4' if typ == 'float' else '>f8'
                nb = n * 3 * np.dtype(dt).itemsize
                out['points'] = np.frombuffer(txt[pos:pos + nb], dtype=dt).reshape(-1, 3).astype(float)
                pos += nb
            else:
                vals = []
                while len(vals) < 3 * n:
                    vals += nextline().split()
                out['points'] = np.array(vals[:3 * n], float).reshape(-1, 3)
        elif key in ('LINES', 'POLYGONS'):
            n, size = int(tok[1]), int(tok[2])
            if len(tok) > 3:
                break
            if binary:
                arr = np.frombuffer(txt[pos:pos + 4 * size], dtype='>i4')
                pos += 4 * size
            else:
                vals = []
                while len(vals) < size:
                    vals += nextline().split()
                arr = np.array(vals[:size], int)
            cells, i = [], 0
            while i < len(arr):
                k = arr[i]
                cells.append(arr[i + 1:i + 1 + k])
                i += 1 + k
            out[key.lower()] = cells
        elif key == 'OFFSETS':      # VTK >= 5.1 legacy: OFFSETS / CONNECTIVITY pairs
            pass
    return out


def _vtk_xml(b):
    """VTK XML PolyData (.vtp): ascii, base64 (optionally zlib) or appended raw/base64 data arrays."""
    root = ET.fromstring(b[:b.find(b'<AppendedData')] + b'</VTKFile>' if b'<AppendedData' in b else b)
    comp = root.get('compressor')
    hdr_t = 'u8' if root.get('header_type') == 'UInt64' else 'u4'
    bo = '<' if root.get('byte_order', 'LittleEndian') == 'LittleEndian' else '>'
    app = None
    if b'<AppendedData' in b:
        s = b.index(b'_', b.index(b'<AppendedData')) + 1
        app = b[s:b.rindex(b'</AppendedData>')]
        app_enc = re.search(rb'<AppendedData[^>]*encoding="(\w+)"', b).group(1).decode()
    tmap = {'Float32': 'f4', 'Float64': 'f8', 'Int32': 'i4', 'Int64': 'i8', 'UInt32': 'u4', 'UInt64': 'u8', 'Int8': 'i1',
            'UInt8': 'u1', 'Int16': 'i2', 'UInt16': 'u2'}

    def decode(da):
        dt = np.dtype(bo + tmap[da.get('type')])
        fmt = da.get('format', 'ascii')
        if fmt == 'ascii':
            return np.array(da.text.split(), dtype=dt)
        hs = np.dtype(bo + hdr_t).itemsize
        if fmt == 'binary':
            raw = base64.b64decode(''.join(da.text.split()))
            data = raw
        else:
            off = int(da.get('offset'))
            data = app[off:] if app_enc == 'raw' else base64.b64decode(app[off:].split(b'_')[0] + b'==')
        if comp:
            nb = int(np.frombuffer(data[:hs], dtype=bo + hdr_t)[0])
            hd = np.frombuffer(data[:hs * (3 + nb)], dtype=bo + hdr_t)
            sizes = hd[3:3 + nb]
            p, chunks = hs * (3 + nb), []
            for s in sizes:
                chunks.append(zlib.decompress(data[p:p + int(s)]))
                p += int(s)
            return np.frombuffer(b''.join(chunks), dtype=dt)
        n = int(np.frombuffer(data[:hs], dtype=bo + hdr_t)[0])
        return np.frombuffer(data[hs:hs + n], dtype=dt)
    piece = root.find('.//Piece')
    out = {}
    pts = piece.find('Points/DataArray')
    out['points'] = decode(pts).astype(float).reshape(-1, 3)
    for tag in ('Lines', 'Polys'):
        el = piece.find(tag)
        if el is None:
            continue
        arrs = {da.get('Name'): decode(da) for da in el.findall('DataArray')}
        if 'connectivity' in arrs and len(arrs['connectivity']):
            con, offs = arrs['connectivity'], arrs['offsets']
            starts = np.r_[0, offs[:-1]]
            out['lines' if tag == 'Lines' else 'polygons'] = [con[s:e] for s, e in zip(starts, offs)]
    return out


def read_poly(name, b):
    n = name.lower()
    if n.endswith('.vtp'):
        return _vtk_xml(b)
    if n.endswith('.vtk'):
        return _vtk_xml(b) if b.lstrip().startswith(b'<') else _vtk_legacy(b)
    if n.endswith('.stl'):
        if b[:5] == b'solid' and b'facet' in b[:400]:
            v = np.array(re.findall(rb'vertex\s+(\S+)\s+(\S+)\s+(\S+)', b), float)
        else:
            nt = struct.unpack('<I', b[80:84])[0]
            rec = np.frombuffer(b[84:84 + nt * 50], dtype=np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]))
            v = rec['v'].reshape(-1, 3).astype(float)
        return {'points': v, 'polygons': [np.arange(i, i + 3) for i in range(0, len(v), 3)]}
    if n.endswith('.mrk.json') or n.endswith('.json'):
        js = json.loads(b)
        lines, pts = [], []
        for mk in js.get('markups', []):
            cp = [c['position'] for c in mk.get('controlPoints', []) if 'position' in c]
            if mk.get('coordinateSystem', 'LPS') == 'LPS':
                cp = [[-x, -y, z] for x, y, z in cp]
            if cp:
                lines.append(np.arange(len(pts), len(pts) + len(cp)))
                pts += cp
        return {'points': np.array(pts, float).reshape(-1, 3), 'lines': lines}
    return None


def lps_flag(b):
    """Slicer writes models in LPS since 4.11 (a '# SPACE=LPS' / 'SPACE=LPS' comment); flip to RAS."""
    head = b[:600]
    return b'SPACE=LPS' in head or b'"LPS"' in head


# --------------------------------------------------------------------------------------------------------------- inspect
def inspect(path):
    z = zipfile.ZipFile(path)
    print(path, len(z.namelist()), 'files')
    for nm in z.namelist():
        info = z.getinfo(nm)
        line = f'  {nm}  ({info.file_size / 1e6:.2f} MB)'
        low = nm.lower()
        try:
            if low.endswith('.mrml'):
                root = ET.fromstring(z.read(nm))
                kinds = {}
                for el in root:
                    kinds.setdefault(el.tag, []).append(el.get('name'))
                print(line)
                for k, v in kinds.items():
                    print(f'      {k}: {len(v)}  {v[:8]}')
                continue
            if low.endswith('.nrrd') or low.endswith('.nhdr'):
                h, _ = read_nrrd_header(z.read(nm)[:20000])
                extra = {k: h[k] for k in ('sizes', 'type', 'kinds', 'space', 'space directions', 'measurement frame') if k in h}
                line += f'  NRRD {extra}'
                segs = sorted({v for k, v in h.items() if re.match(r'segment\d+_name$', k)})
                if segs:
                    line += f'  segments: {segs}'
            else:
                d = read_poly(nm, z.read(nm))
                if d:
                    line += (f'  points {len(d["points"])}, lines {len(d.get("lines", []))}, '
                             f'polygons {len(d.get("polygons", []))}, LPS {lps_flag(z.read(nm))}')
                    if len(d['points']):
                        P = d['points']
                        line += f'  bbox {P.min(0).round(1)} .. {P.max(0).round(1)}'
        except Exception as ex:  # noqa: BLE001
            line += f'  (not read: {ex})'
        print(line)


# ------------------------------------------------------------------------------------------------------------------- map
def polylines(z, pattern):
    P, D = [], []
    for nm in z.namelist():
        if not re.search(pattern, nm, re.I):
            continue
        b = z.read(nm)
        d = read_poly(nm, b)
        if not d or not d.get('lines'):
            continue
        X = d['points'] * (np.array([-1, -1, 1]) if lps_flag(b) else 1)
        for ln in d['lines']:
            if len(ln) < 2:
                continue
            a, c = X[ln[:-1]], X[ln[1:]]
            v = c - a
            L = np.linalg.norm(v, axis=1)
            ok = L > 1e-9
            P.append((a + c)[ok] / 2)
            D.append(v[ok] / L[ok, None])
        print(f'fibres: {nm}: {len(d["lines"])} lines')
    return np.vstack(P), np.vstack(D)


def surface_points(z, pattern):
    pts = []
    for nm in z.namelist():
        if re.search(pattern, nm, re.I):
            b = z.read(nm)
            d = read_poly(nm, b)
            if d is not None and len(d['points']):
                pts.append(d['points'] * (np.array([-1, -1, 1]) if lps_flag(b) else 1))
                print(f'levator surface: {nm}: {len(d["points"])} points')
    return np.vstack(pts)


def rigid_icp(src, dst, iters=60):
    """Rigid ICP src -> dst from 8 principal-axis starts; returns R, t, rms (src @ R.T + t ~ dst)."""
    cs, cd = src.mean(0), dst.mean(0)
    _, _, Vs = np.linalg.svd(src - cs, full_matrices=False)
    _, _, Vd = np.linalg.svd(dst - cd, full_matrices=False)
    T = cKDTree(dst)
    best = None
    for sx in (1, -1):
        for sy in (1, -1):
            S = np.diag([sx, sy, sx * sy * np.linalg.det(Vs) * np.linalg.det(Vd)])
            R = Vd.T @ S @ Vs
            t = cd - R @ cs
            for _ in range(iters):
                Y = src @ R.T + t
                d, j = T.query(Y)
                keep = d <= np.percentile(d, 90)
                A, B = src[keep], dst[j[keep]]
                ca, cb = A.mean(0), B.mean(0)
                U, _, Vt = np.linalg.svd((A - ca).T @ (B - cb))
                Dm = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
                R = Vt.T @ Dm @ U.T
                t = cb - R @ ca
            rms = float(np.sqrt(np.mean(T.query(src @ R.T + t)[0] ** 2)))
            if best is None or rms < best[2]:
                best = (R, t, rms)
    return best


def la_elements(feb):
    r = ET.parse(feb).getroot()
    m = r.find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in m.findall('Nodes') for n in b}
    E = {e.get('name'): [[int(v) for v in x.text.split(',')] for x in e] for e in m.findall('Elements')
         if e.get('name') in LA_DOMAINS}
    return X, E


def fibre_field(cent, nrm, P, D, max_dist):
    T = cKDTree(P)
    out = np.zeros((len(cent), 3))
    coh = np.zeros(len(cent))
    cnt = np.zeros(len(cent), int)
    for i, (c, n) in enumerate(zip(cent, nrm)):
        idx = T.query_ball_point(c, max_dist)
        if not idx:
            continue
        d = np.linalg.norm(P[idx] - c, axis=1)
        w = 1 / np.maximum(d, 0.5)
        v = D[idx] - np.outer(D[idx] @ n, n)          # into the tangent plane
        L = np.linalg.norm(v, axis=1)
        ok = L > 1e-6
        if not ok.any():
            continue
        v = v[ok] / L[ok, None]
        M = (w[ok, None, None] * v[:, :, None] * v[:, None, :]).sum(0) / w[ok].sum()
        ev, evec = np.linalg.eigh(M)
        out[i] = evec[:, -1]
        coh[i] = ev[-1] - ev[-2]
        cnt[i] = ok.sum()
    miss = cnt == 0
    if miss.any() and (~miss).any():
        Tc = cKDTree(cent[~miss])
        _, j = Tc.query(cent[miss])
        src = np.where(~miss)[0][j]
        for i, s in zip(np.where(miss)[0], src):
            v = out[s] - (out[s] @ nrm[i]) * nrm[i]
            out[i] = v / np.linalg.norm(v)
    return out, coh, cnt


def run_map(mrb, feb, fib_pat, la_pat, out, max_dist):
    z = zipfile.ZipFile(mrb)
    P, D = polylines(z, fib_pat)
    S = surface_points(z, la_pat)
    X, E = la_elements(feb)
    la_nodes = np.array([X[n] for e in E.values() for el in e for n in el])
    scale = np.sqrt(((S - S.mean(0)) ** 2).sum(1).mean()) / np.sqrt(((la_nodes - la_nodes.mean(0)) ** 2).sum(1).mean())
    R, t, rms = rigid_icp(S, np.unique(la_nodes, axis=0))
    print(f'registration: rigid ICP rms {rms:.2f} mm (scene / model size ratio {scale:.2f}: ~1 expected)')
    P = P @ R.T + t
    D = D @ R.T
    blocks, rows = [], ['element_lid,domain,cx,cy,cz,ax,ay,az,coherence,n_segments']
    for dom, els in E.items():
        C, N = [], []
        for el in els:
            Q = np.array([X[n] for n in el])
            v = np.cross(Q[2] - Q[0], Q[-1] - Q[1]) if len(el) == 4 else np.cross(Q[1] - Q[0], Q[2] - Q[0])
            C.append(Q.mean(0))
            N.append(v / np.linalg.norm(v))
        C, N = np.array(C), np.array(N)
        A, coh, cnt = fibre_field(C, N, P, D, max_dist)
        Dn = np.cross(N, A)
        lines = [f'\t\t<ElementData type="mat_axis" elem_set="{dom}">']
        for k in range(len(els)):
            lines.append(f'\t\t\t<elem lid="{k + 1}"><a>{A[k][0]:.6f},{A[k][1]:.6f},{A[k][2]:.6f}</a>'
                         f'<d>{Dn[k][0]:.6f},{Dn[k][1]:.6f},{Dn[k][2]:.6f}</d></elem>')
            rows.append(f'{k + 1},{dom},{C[k][0]:.3f},{C[k][1]:.3f},{C[k][2]:.3f},{A[k][0]:.5f},{A[k][1]:.5f},'
                        f'{A[k][2]:.5f},{coh[k]:.3f},{cnt[k]}')
        lines.append('\t\t</ElementData>')
        blocks.append('\n'.join(lines))
        print(f'{dom}: {len(els)} elements, {np.mean(cnt > 0) * 100:.0f} % with fibres within {max_dist} mm, '
              f'median coherence {np.median(coh[cnt > 0]) if (cnt > 0).any() else 0:.2f}')
    open(out, 'w').write('\n'.join(blocks) + '\n')
    open(os.path.splitext(out)[0] + '.csv', 'w').write('\n'.join(rows) + '\n')
    print('wrote', out, 'and the .csv')


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('inspect', 'map'):
        sys.exit(__doc__)
    if sys.argv[1] == 'inspect':
        inspect(sys.argv[2])
    else:
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument('cmd'); ap.add_argument('mrb'); ap.add_argument('run')
        ap.add_argument('--fibres', required=True); ap.add_argument('--la-surface', required=True)
        ap.add_argument('--out', default=None); ap.add_argument('--max-dist', type=float, default=5)
        a = ap.parse_args()
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from paths import RUNS_DIR
        feb = a.run if a.run.endswith('.feb') else os.path.join(RUNS_DIR, a.run, a.run + '.feb')
        run_map(a.mrb, feb, a.fibres, a.la_surface, a.out or os.path.splitext(a.mrb)[0] + '_LA_mat_axis.xml', a.max_dist)
