"""Batch 91 (2026-09-27 night): the PVW edge on the LA with the faithful walls is a USL pinch (NOT IN SOURCE diagnostic).

In L85_newline_vwyeoh_rhoi0 (t 0.70-0.73) the PVW nodes that reverse direction on the LA step to step (0.2-0.8 m/s,
0.005-0.015 MPa) are 1892-1898 (left) and 1984-1988 (right): uterosacral connector ends (USL-L_conn_2, USL-R_conn_2) on
the lateral edge row of the PVW_LA primary surface, the same mechanism as the P-arcus pinch the user's option (b) removed
(the connectors pull the PVW's edge onto the LA at a point). The USL families hold 22 nodes of that edge (USL-L/R_conn_1:
3 each; USL-L/R_conn_2: 8 each). One change from L85_newline_vwyeoh_rhoi0:
  L86_newline_vwyeoh_rhoi0_uslfacets   the PVW_LA primary facets holding a USL connector end removed (NOT IN SOURCE)
usage: py -3.10 build_batch91.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WANT = set(sys.argv[1:])


def facets_out(m, prefix='USL-', surface='PVW_LA_primary'):
    """NOT IN SOURCE: drop every <surface> facet holding an end of a DiscreteSet whose name starts with prefix."""
    mesh = m.mesh
    sets = [d.get('name') for d in mesh.findall('DiscreteSet') if d.get('name').startswith(prefix)]
    pn = {int(v) for d in mesh.findall('DiscreteSet') if d.get('name') in sets for e in d for v in e.text.split(',')}
    s = next(s for s in mesh.findall('Surface') if s.get('name') == surface)
    n0 = len(s)
    gone = [f for f in list(s) if {int(v) for v in f.text.split(',')} & pn]
    for f in gone:
        s.remove(f)
    for i, f in enumerate(s, 1):
        f.set('id', str(i))
    m.log.append(f'NOT IN SOURCE (PVW-on-LA pinch diagnostic): the {len(gone)} {surface} facets holding an end of '
                 f'{", ".join(sets)} removed from the contact surface ({len(s)} of {n0} left), as option (b) did for P-arcus')
    return len(gone)


BUILDS = (('L86_newline_vwyeoh_rhoi0_uslfacets', 'L85_newline_vwyeoh_rhoi0', facets_out),)
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(os.path.join(RUNS, base, base + '.feb'))
        fn(m)
        print(m.log[-1])
        emit(name, m)
