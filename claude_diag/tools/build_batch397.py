"""Batch 397 (2026-10-22 morning; the user's two ideas against the arch, the canal rising instead of prolapsing out; 'not for
the long term, but I'd like to see what happens'). Both NOT IN SOURCE, test devices only.
  --plane : a fixed rigid plane above the apex (normal -z, facing the model) at GAP mm above the highest rest z of the solid
            domains, with frictionless sliding-elastic contact (PENALTY PEN auto, two_pass 0, search_radius SR) against the
            outward boundary faces of all solid domains whose outward normal has n_z > 0.3 and centroid within ZBAND mm below
            the plane: the parts cannot move cranially (+z) past it, only slide forward / back / sideways (the user: "so it
            physically can't go up. It can only slide back or slide forwards"). Syntax checked in claude_diag/mini/rigidplane
            (penalty 100 held a block at the plane, it slid freely).
  --push  : a downward traction (0, 0, -1) x P MPa on Load-AVW_surf (the AVW wall-pressure faces) on its own load curve
            0 at t 0 -> 1 at T1 -> 0 at T2 (SMOOTH), then 0 (the user: "a downward force on the avw that decreases and
            then goes to 0 late in the run"). Syntax checked in claude_diag/mini/traction.
usage: py -3.10 build_batch397.py SRC_RUN NEW_RUN --plane [--gap 0.5] [--pen 30] [--sr 1] [--zband 12]
       py -3.10 build_batch397.py SRC_RUN NEW_RUN --push [--p 0.002] [--t1 0.5] [--t2 0.9]
       py -3.10 build_batch397.py SRC_RUN NEW_RUN --press [--p 0.002] [--t1 0.5] [--t2 0.8]   (L401: extra AVW pressure)"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402
from build_batch297 import HF  # noqa: E402
from junction_view_data import load  # noqa: E402

WF = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]   # penta6 outward faces


def add_load_curve(t, name, pts):
    ids = [int(v) for v in re.findall(r'<load_controller id="(\d+)"', t)]
    lc = max(ids) + 1
    blk = (f'\t\t<load_controller id="{lc}" name="{name}" type="loadcurve">\n\t\t\t<interpolate>SMOOTH</interpolate>\n'
           '\t\t\t<extend>CONSTANT</extend>\n\t\t\t<points>\n' + ''.join(f'\t\t\t\t<pt>{a:g},{b:g}</pt>\n' for a, b in pts) +
           '\t\t\t</points>\n\t\t</load_controller>\n')
    i = t.index('\t</LoadData>')
    return t[:i] + blk + t[i:], lc


def plane(src, new, gap, pen, sr, zband):
    t, X, E = load(src)
    solid = {k: v for k, v in E.items() if len(v[1]) in (6, 8)}
    cnt, own = {}, {}
    for nm, c in solid.values():
        for f in (HF if len(c) == 8 else WF):
            q = tuple(c[i] for i in f); k = frozenset(q)
            cnt[k] = cnt.get(k, 0) + 1; own[k] = (q, nm)
    faces = [own[k] for k, v in cnt.items() if v == 1]
    ztop = max(X[v][2] for _, c in solid.values() for v in c)
    zp = ztop + gap
    prim, doms = [], {}
    for q, nm in faces:
        P = np.array([X[v] for v in q]); c = P.mean(0)
        n = np.cross(P[1] - P[0], P[2] - P[0]) if len(q) == 3 else np.cross(P[2] - P[0], P[3] - P[1])
        n /= np.linalg.norm(n)
        if n[2] > 0.3 and c[2] > zp - zband:
            prim.append(q); doms[nm] = doms.get(nm, 0) + 1
    allP = np.array([X[v] for q in prim for v in q])
    x0, x1 = allP[:, 0].min() - 40, allP[:, 0].max() + 40
    y0, y1 = allP[:, 1].min() - 60, allP[:, 1].max() + 120
    nid = max(X) + 1
    eid = max(E) + 1
    nx = ny = 4
    xs, ys = np.linspace(x0, x1, nx + 1), np.linspace(y0, y1, ny + 1)
    pid = {}
    nodes = '\t\t<Nodes name="cap_plane_nodes">\n'
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            pid[i, j] = nid; nodes += f'\t\t\t<node id="{nid}">{x:.4f},{y:.4f},{zp:.4f}</node>\n'; nid += 1
    nodes += '\t\t</Nodes>\n'
    # winding (i,j) -> (i,j+1) -> (i+1,j+1) -> (i+1,j): normal -z (faces the model below), checked in the mini test
    quads = [(pid[i, j], pid[i, j + 1], pid[i + 1, j + 1], pid[i + 1, j]) for j in range(ny) for i in range(nx)]
    els = '\t\t<Elements type="quad4" name="cap_plane">\n' + ''.join(
        f'\t\t\t<elem id="{eid + k}">{",".join(map(str, q))}</elem>\n' for k, q in enumerate(quads)) + '\t\t</Elements>\n'
    i = t.rindex('\t\t</Nodes>\n') + len('\t\t</Nodes>\n')
    t = t[:i] + nodes + t[i:]
    i = t.rindex('\t\t</Elements>\n') + len('\t\t</Elements>\n')
    t = t[:i] + els + t[i:]
    surf = ('\t\t<Surface name="cap_plane_primary">\n' + ''.join(
        f'\t\t\t<{"tri3" if len(q) == 3 else "quad4"} id="{k}">{",".join(map(str, q))}</{"tri3" if len(q) == 3 else "quad4"}>\n'
        for k, q in enumerate(prim, 1)) + '\t\t</Surface>\n\t\t<Surface name="cap_plane_secondary">\n' + ''.join(
        f'\t\t\t<quad4 id="{k}">{",".join(map(str, q))}</quad4>\n' for k, q in enumerate(quads, 1)) +
        '\t\t</Surface>\n\t\t<SurfacePair name="cap_plane">\n\t\t\t<primary>cap_plane_primary</primary>\n'
        '\t\t\t<secondary>cap_plane_secondary</secondary>\n\t\t</SurfacePair>\n')
    i = t.index('\t\t<SurfacePair name=')
    t = t[:i] + surf + t[i:]
    mids = [int(v) for v in re.findall(r'<material id="(\d+)"', t)]
    mat = (f'\t\t<material id="{max(mids) + 1}" name="cap_plane_rigid" type="rigid body">\n'
           '\t\t\t<density>1.06e-09</density>\n\t\t</material>\n')
    i = t.index('</Material>')
    t = t[:i] + mat + t[i:]
    i = t.index('\t</MeshDomains>')
    t = t[:i] + '\t\t<ShellDomain name="cap_plane" mat="cap_plane_rigid">\n\t\t\t<shell_thickness>0.1</shell_thickness>\n\t\t</ShellDomain>\n' + t[i:]
    assert '<Rigid>' not in t
    rig = ('\t<Rigid>\n\t\t<rigid_bc name="cap_plane_fixed" type="rigid_fixed">\n\t\t\t<rb>cap_plane_rigid</rb>\n' +
           ''.join(f'\t\t\t<{d}_dof>1</{d}_dof>\n' for d in ('Rx', 'Ry', 'Rz', 'Ru', 'Rv', 'Rw')) + '\t\t</rigid_bc>\n\t</Rigid>\n')
    i = t.index('\t</Boundary>\n') + len('\t</Boundary>\n')
    t = t[:i] + rig + t[i:]
    a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
    c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', 'name="cap_plane" surface_pair="cap_plane"')
         .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
         .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
    for s_ in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>', '<fric_coeff>0</fric_coeff>'):
        assert s_ in c, s_
    t = t[:b] + '\t\t' + c + t[b:]
    msg = (f'NOT IN SOURCE (test device, the user 10-22): a fixed rigid plane (cap_plane, rigid body fixed in 6 DOFs, 4 x 4 quad4, '
           f'normal -z) at z {zp:.2f} = the highest solid rest z {ztop:.2f} + {gap:g} mm, x {x0:.0f}..{x1:.0f}, y {y0:.0f}..{y1:.0f}; '
           f'frictionless sliding-elastic cap_plane (PENALTY {pen:g} auto, two_pass 0, search_radius {sr:g}, offset 0, as walls_LA '
           f'otherwise) against {len(prim)} outward boundary faces (n_z > 0.3, centroid within {zband:g} mm below the plane): '
           + ', '.join(f'{k} {v}' for k, v in sorted(doms.items())) + '; the parts cannot move cranially past the plane, only slide')
    return t, msg


def press(src, new, p, t1, t2):
    """(L401, the user 10-22: 'a force on the surface of the avw just like the normal force ... goes away later in the run,
    around 0.8') an extra PRESSURE (normal to the faces, as Load-AVW) on Load-AVW_surf on its own load curve."""
    t = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', src, src + '.feb'), encoding='utf-8').read()
    assert '<Surface name="Load-AVW_surf">' in t and 'name="Load-AVW" surface="Load-AVW_surf" type="pressure"' in t
    t, lc = add_load_curve(t, 'avw_extra', [(0, 0), (t1, 1), (t2, 0)])
    blk = (f'\t\t<surface_load name="avw_extra" surface="Load-AVW_surf" type="pressure">\n\t\t\t<pressure lc="{lc}">{p:g}</pressure>\n'
           '\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n\t\t\t<linear>0</linear>\n\t\t\t<shell_bottom>0</shell_bottom>\n\t\t</surface_load>\n')
    i = t.rindex('</surface_load>\n') + len('</surface_load>\n')
    t = t[:i] + blk + t[i:]
    msg = (f'NOT IN SOURCE (test device, the user 10-22): avw_extra = an extra pressure {p:g} MPa on Load-AVW_surf (as Load-AVW, '
           f'0.014; x{p / 0.014:.2f}; net at rest ~(0, -24.9, -16.9) N x {p / 0.014:.2f}), follower, on load curve {lc} (SMOOTH 0 at '
           f't 0 -> 1 at {t1:g} -> 0 at {t2:g}, then 0)')
    return t, msg


def push(src, new, p, t1, t2):
    t = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', src, src + '.feb'), encoding='utf-8').read()
    assert '<Surface name="Load-AVW_surf">' in t
    t, lc = add_load_curve(t, 'avw_push', [(0, 0), (t1, 1), (t2, 0)])
    blk = (f'\t\t<surface_load name="avw_push" surface="Load-AVW_surf" type="traction">\n\t\t\t<scale lc="{lc}">{p:g}</scale>\n'
           '\t\t\t<traction>0,0,-1</traction>\n\t\t</surface_load>\n')
    i = t.rindex('</surface_load>\n') + len('</surface_load>\n')
    t = t[:i] + blk + t[i:]
    a = t.index('<Surface name="Load-AVW_surf">'); b = t.index('</Surface>', a)
    n = len(re.findall(r'<quad4|<tri3', t[a:b]))
    msg = (f'NOT IN SOURCE (test device, the user 10-22): avw_push = a downward traction (0, 0, -1) x {p:g} MPa on Load-AVW_surf '
           f'({n} faces, ~2299 mm2 -> ~{p * 2299:.1f} N at the peak) on load curve {lc} (SMOOTH 0 at t 0 -> 1 at {t1:g} -> 0 at '
           f'{t2:g}, then 0): holds the AVW / apex down while the arch would form, gone before t 1')
    return t, msg


if __name__ == '__main__':
    a = sys.argv[1:]
    src, new = a[0], a[1]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    if '--press' in a:
        t, msg = press(src, new, g('--p', 0.002), g('--t1', 0.5), g('--t2', 0.8))
    elif '--plane' in a:
        t, msg = plane(src, new, g('--gap', 0.5), g('--pen', 30), g('--sr', 1), g('--zband', 12))
    else:
        t, msg = push(src, new, g('--p', 0.002), g('--t1', 0.5), g('--t2', 0.9))
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch397.py)\n{msg}\n')
