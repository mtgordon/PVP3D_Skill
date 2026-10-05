"""The PVW_LA contact onset in an existing run (Test B of NEXT_SESSION_PROMPT_2026-09-26b.md): for every converged
state in a time window, the facets of the pair that carry pressure (which, how many, first-touch time), their gap, the
speeds of the wall nodes, and the solver's effort per step from the log (step size, iterations, failed attempts).

The plot file's surface section (not parsed by xplt_reader.py) gives each plotted surface's name and facets, so the
surface data can be tied to PVW_LA_primary / PVW_LA_secondary by name.
usage: py -3.10 pvw_la_onset.py RUN T0 T1 [--every K]
"""
import argparse
import os
import re
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'skill', 'abaqus-febio-fea-pipeline', 'scripts'))
from xplt_reader import Xplt, chunks, _name  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
PLT_SURFACE_SECTION = 0x01043000
PLT_SURFACE = 0x01043100
PLT_SURFACE_HDR = 0x01043101
PLT_SURFACE_ID = 0x01043102
PLT_SURFACE_FACES = 0x01043103
PLT_SURFACE_NAME = 0x01043104
PLT_SURFACE_MAX_FACET_NODES = 0x01043105
PLT_FACE_LIST = 0x01043200
PLT_FACE = 0x01043201


class XpltS(Xplt):
    """Xplt plus the mesh's surface section: self.surfaces = [{'id', 'name', 'nf', 'faces': [(fid, [node idx...])]}]."""

    def _read_mesh(self, buf):
        super()._read_mesh(buf)
        self.surfaces = []
        for cid, o, s in chunks(buf):
            if cid != PLT_SURFACE_SECTION:
                continue
            for c2, o2, s2 in chunks(buf, o, o + s):
                if c2 != PLT_SURFACE:
                    continue
                srf = {'faces': []}
                for c3, o3, s3 in chunks(buf, o2, o2 + s2):
                    if c3 == PLT_SURFACE_HDR:
                        for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                            if c4 == PLT_SURFACE_ID:
                                srf['id'] = struct.unpack_from('<I', buf, o4)[0]
                            elif c4 == PLT_SURFACE_FACES:
                                srf['nf'] = struct.unpack_from('<I', buf, o4)[0]
                            elif c4 == PLT_SURFACE_NAME:
                                srf['name'] = _name(buf[o4:o4 + s4])
                            elif c4 == PLT_SURFACE_MAX_FACET_NODES:
                                srf['maxn'] = struct.unpack_from('<I', buf, o4)[0]
                    elif c3 == PLT_FACE_LIST:
                        for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                            if c4 == PLT_FACE:
                                a = np.frombuffer(buf, dtype=np.int32, count=s4 // 4, offset=o4)
                                srf['faces'].append((int(a[0]), [int(v) for v in a[2:2 + int(a[1])]]))
                self.surfaces.append(srf)


def log_steps(path):
    """[(t, dt, iterations, failed attempts before it)] per converged step, from the log."""
    txt = open(path, encoding='latin-1').read()
    out, prev_t, fails, its = [], 0.0, 0, 0
    for m in re.finditer(r'^------- (converged|failed to converge) at time : ([0-9.eE+-]+)|'
                         r'^\s+number of iterations\s+:\s+(\d+)', txt, re.M):
        if m.group(3):
            its = int(m.group(3))
        elif m.group(1) == 'converged':
            t = float(m.group(2))
            out.append((t, t - prev_t, its, fails))
            prev_t, fails = t, 0
        else:
            fails += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('t0', type=float); ap.add_argument('t1', type=float)
    ap.add_argument('--every', type=int, default=1)
    a = ap.parse_args()
    x = XpltS(os.path.join(RUNS, a.run, a.run + '.xplt'))
    names = [s.get('name') for s in x.surfaces]
    print('plotted surfaces:', [(s.get('id'), s.get('name'), s.get('nf')) for s in x.surfaces])
    print('surface variables:', [(it.get('name'), it.get('fmt')) for it in x.dict['surface']])
    sp = {k: next(s for s in x.surfaces if s.get('name') == 'PVW_LA_' + k) for k in ('primary', 'secondary')}
    # converged states: status flag 0 if present, else by the log's converged times
    conv = [i for i, s in enumerate(x.states) if s[1] == 0] or list(range(len(x.states)))
    steps = log_steps(os.path.join(RUNS, a.run, a.run + '.log'))
    st_t = np.array([s[0] for s in steps])
    win = [k for k in conv if a.t0 <= x.states[k][0] <= a.t1][::a.every]
    first_touch = {'primary': {}, 'secondary': {}}
    pvw_nodes = sorted({n for _, f in sp['primary']['faces'] for n in f})
    la_nodes = sorted({n for _, f in sp['secondary']['faces'] for n in f})
    print(f"primary {sp['primary']['nf']} facets / {len(pvw_nodes)} nodes; secondary {sp['secondary']['nf']} facets / "
          f"{len(la_nodes)} nodes")
    hdr = ('t', 'dt', 'its', 'fails', 'prim n>0', 'prim pmax', 'sec n>0', 'sec pmax', 'prim gap min', 'PVW vmax mm/s',
           'LA vmax mm/s', 'touching (primary facet ids)')
    print(' | '.join(hdr))
    for k in win:
        t = x.states[k][0]
        u = x.var(k, 'displacement')
        # node speed over the step, from the displacement change since the previous converged state (the plotted
        # "velocity" is an element variable in these files)
        kp = conv[conv.index(k) - 1] if conv.index(k) > 0 else None
        up, tp = (x.var(kp, 'displacement'), x.states[kp][0]) if kp is not None else (None, None)
        v = (u - up) / (t - tp) if up is not None and t > tp else None
        j = int(np.argmin(np.abs(st_t - t))) if len(st_t) else None
        dt, its, fails = (steps[j][1], steps[j][2], steps[j][3]) if j is not None and abs(st_t[j] - t) < 1e-5 else ('', '', '')
        cp = x.var(k, 'contact pressure') or {}
        gp = x.var(k, 'contact gap') or {}
        row = [f'{t:.5f}', f'{dt:.2e}' if dt != '' else '', str(its), str(fails)]
        touching = {}
        for side in ('primary', 'secondary'):
            sid = sp[side]['id']
            p = np.asarray(cp.get(sid, [])).reshape(sp[side]['nf'], -1).max(axis=1) if sid in cp else np.zeros(0)
            on = np.where(p > 1e-9)[0]
            touching[side] = on
            for f in on:
                first_touch[side].setdefault(int(f), t)
            row += [str(len(on)), f'{p.max():.3g}' if p.size else '']
        g = np.asarray(gp.get(sp['primary']['id'], [])).reshape(sp['primary']['nf'], -1) if sp['primary']['id'] in gp else None
        row.append(f'{g.min():.3g}' if g is not None and g.size else '')
        sp_v = np.linalg.norm(v, axis=1) if v is not None else None
        row.append(f'{sp_v[pvw_nodes].max():.0f}' if sp_v is not None else '')
        row.append(f'{sp_v[la_nodes].max():.0f}' if sp_v is not None else '')
        row.append(','.join(str(sp['primary']['faces'][f][0]) for f in touching['primary'][:8]))
        print(' | '.join(row))
    for side in ('primary', 'secondary'):
        ft = sorted(first_touch[side].items(), key=lambda kv: kv[1])
        print(f'{side}: {len(ft)} facets ever touching in the window; first touches:',
              [(sp[side]['faces'][f][0], [int(x.node_ids[n]) for n in sp[side]['faces'][f][1]], round(t, 5))
               for f, t in ft[:10]])


if __name__ == '__main__':
    main()
