"""Are any FEBio shell elements pinched or inverted *before any load is applied*?

FEBio shells are solid-like: each element integrates through its thickness along
nodal directors, which by default (shell_normal_nodal = 1) are the average of the
normals of the shell elements around each node. Where a mesh folds more tightly
than the shell thickness allows, the offset faces of the elements at the fold
cross over, and the element starts out (nearly) inverted. It then fails at a tiny
load fraction with "negative jacobian". Tuning materials, penalties or solver
settings cannot fix that. Abaqus conventional shells (S4R/S3R) never form this
through-thickness volume, so a thick shell imported from Abaqus can carry the
defect silently.

**Geometry convention (this matters -- an earlier version of this script got it
wrong and roughly halved the apparent pinching on real models):** FEBio 4's
default shell places the mesh NODES at the shell's TOP face; the element extends
a full thickness t along the *negative* averaged nodal normal to reach the
bottom face (FEBio User Manual 3.6.2.2). Only `<shell_formulation>0</shell_formulation>`
in `<Control>` (the pre-2.6 legacy mode, compatible-strain formulations only)
recovers the older mid-surface convention (nodes at the geometric middle,
+-t/2 either side) -- inferred from the manual's prose, not verified against a
live example, so treat the auto-detect as a best-effort default and pass
--mid-surface to force it if a model is known to use shell_formulation=0. This
script auto-detects `<shell_formulation>` and otherwise assumes the FEBio 4
top-face default. Getting this backwards makes the "top" side trivially
uninvertible (ratio always 1.000, since that's exactly where the reference
nodes are) and roughly halves the offset that actually matters on the other
side: on one real 4 mm shell, the true Gauss-point-inverted count was 86
elements, and the old (mid-surface) offset formula found only 1.

For every selected shell element this computes, along the averaged nodal
normals, the signed area of the element's offset surface at the two faces and
at the two through-thickness Gauss points, divided by the reference
(node-surface) area. 1 = prism-like, < 0.25 = badly pinched, <= 0 = inverted.
Under the FEBio-4 default convention the "top" face IS the node surface and is
trivially always 1.000 -- only "bottom" and the Gauss points can show
pinching; under --mid-surface both faces are informative. A domain whose
<ShellDomain> already has shell_normal_nodal=0 is scored against its OWN
per-element normal, not a neighbor-averaged one (matching what FEBio actually
does for it) -- such a domain will always come back ratio=1.000 everywhere,
which is the expected way to confirm the fix took, not a sign the check did
nothing. This is a geometric *screen*, not a copy of FEBio's Jacobian: it
assumes flat facets and unit averaged normals, so treat flagged elements as
suspects, not a certainty. In
the reference session the worst-ranked elements (face ratio well below 0 on a
4 mm shell) were exactly the ones that inverted in every failing run, the
model still started normally, and <shell_normal_nodal>0</shell_normal_nodal>
on those ShellDomains fixed it (element normals keep each element a proper
prism; thickness unchanged). Confirm with a run, and compare --per-domain vs
the default when domains meet at a fold, since the two averaging scopes can
rank elements differently.

usage: py -3 shell_pinch_check.py model.feb [--domains A,B,...] [--thickness T] [--per-domain] [--mid-surface]
  --domains      comma-separated ShellDomain names (default: every ShellDomain)
  --thickness    override the <shell_thickness> read from each ShellDomain
  --per-domain   average nodal normals within each domain only (default: across all
                 selected domains). Restricting --domains to the thick shell(s) of
                 interest gives the cleanest ranking; thin shells mostly show up only
                 where they have genuine slivers.
  --mid-surface  force the pre-2.6 mid-surface convention instead of auto-detecting
                 <shell_formulation> from <Control> (FEBio 4 default: top-face).
Needs Python 3 + numpy; imports feb_model.py from this directory.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feb_model import Feb  # noqa: E402


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
    ap.add_argument('--thickness', type=float)
    ap.add_argument('--per-domain', action='store_true')
    ap.add_argument('--mid-surface', action='store_true',
                     help='force the pre-2.6 mid-surface convention instead of auto-detecting '
                          '<shell_formulation> from <Control>')
    a = ap.parse_args()

    f = Feb(a.feb)
    ctrl = f.section('Control')
    sf = ctrl.find('shell_formulation') if ctrl is not None else None
    mid_surface = a.mid_surface or (sf is not None and sf.text.strip() == '0')

    thick, nodal_flag = {}, {}
    for dom in f.section('MeshDomains'):
        if dom.tag == 'ShellDomain':
            t = dom.find('shell_thickness')
            thick[dom.get('name')] = a.thickness or (float(t.text.split(',')[0]) if t is not None else 0.0)
            sn = dom.find('shell_normal_nodal')
            nodal_flag[dom.get('name')] = (sn is None) or sn.text.strip() != '0'
    names = a.domains.split(',') if a.domains else sorted(thick)
    names = [n for n in names if n in f.elem_blocks and n in thick]
    if not names:
        sys.exit('no ShellDomains selected')

    elems = []
    for n in names:
        for eid, conn in f.elem_blocks[n][1].items():
            P = np.array([f.nodes[i] for i in conn])
            elems.append((n, eid, conn, P, elem_normal(P)))
    acc = defaultdict(lambda: np.zeros(3))
    for n, eid, conn, P, nn in elems:
        if nodal_flag.get(n, True):    # only accumulate where the domain actually nodal-averages
            for i in conn:
                acc[(n if a.per_domain else '*', i)] += nn

    # Face/Gauss-point offsets, in units of thickness t, measured from the node along -D:
    #   mid-surface: nodes at the geometric mid-surface, offsets +-t/2 either side.
    #   top-face (FEBio 4 default): nodes ARE the top face (offset 0, trivially ratio 1),
    #     the shell extends a full t along -D to the bottom face.
    if mid_surface:
        s_face = (0.5, -0.5)
        s_gp = (0.5 / np.sqrt(3), -0.5 / np.sqrt(3))
    else:
        s_face = (0.0, 1.0)
        s_gp = ((1 - 1 / np.sqrt(3)) / 2, (1 + 1 / np.sqrt(3)) / 2)

    rows = []
    for n, eid, conn, P, nn in elems:
        t = thick[n]
        if nodal_flag.get(n, True):
            D = np.array([acc[(n if a.per_domain else '*', i)] for i in conn])
            D /= np.linalg.norm(D, axis=1)[:, None]
        else:
            # shell_normal_nodal=0: FEBio uses the element's OWN normal, not a neighbor average --
            # a single flat element offset along its own normal can never appear pinched (ratio
            # stays exactly 1), which is the point of the fix this screen is meant to confirm.
            D = np.tile(nn, (len(conn), 1))
        a0 = signed_area(P, nn)
        faces = [signed_area(P - s * t * D, nn) / a0 for s in s_face]
        gps = [signed_area(P - s * t * D, nn) / a0 for s in s_gp]
        edges = [np.linalg.norm(P[k] - P[(k + 1) % len(P)]) for k in range(len(P))]
        rows.append((min(faces), n, eid, faces, min(edges), t, min(gps)))

    r = np.array([x[0] for x in rows])
    g = np.array([x[6] for x in rows])
    print(f'{len(rows)} shell elements in {names}')
    print(f'  convention: {"mid-surface (+-t/2)" if mid_surface else "top-face (FEBio 4 default: nodes = top, extends -t)"}')
    print(f'  shell_normal_nodal off (element normals) in: '
          f'{[n for n in names if not nodal_flag.get(n, True)] or "none"}')
    print(f'  face inverted at rest (face ratio <= 0): {int((r <= 0).sum())}   pinched < 0.25: {int((r < 0.25).sum())}'
          f'   < 0.5: {int((r < 0.5).sum())}   min face ratio {r.min():.3f}')
    print(f'  Gauss-point level also inverted (ratio <= 0): {int((g <= 0).sum())}   min GP ratio {g.min():.3f}')
    print('  worst elements (domain, id, face ratios top/bottom, min edge, thickness):')
    for x in sorted(rows)[:15]:
        print(f'    {x[1]:24s} {x[2]:8d}  faces {x[3][0]:7.3f} {x[3][1]:7.3f}   min edge {x[4]:.2f}   t {x[5]:g}')
    if (r <= 0).any():
        print('  -> suspects: try <shell_normal_nodal>0</shell_normal_nodal> on the affected thick-shell '
              'domains and re-run (domains that already have it are not pinched by averaging)')


if __name__ == '__main__':
    main()
