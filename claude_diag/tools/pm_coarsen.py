"""Remesh the PM shell across its width (2026-10-08 night, unattended; the user: "speed up the runs ... you can also try
remeshing"). The slow low-Parcus runs fail mostly on negative Jacobians in PM_Plane (L226_cl30p062_diag_P1: 38 of the first 42
reports), in a cluster of shell quads only 0.36-0.68 mm wide across the band (aspect 7-8) under a 2 mm shell thickness.
The PM (NOT IN SOURCE: the user's PM structure, 2026-10-02) is a flat structured band: 11 rows x 36 nodes (ids base + 36 r + c),
35 x 10 quads; row 0 = the clamped outer arc (BC-PM-outer-arc), rows 9-10 carry the PM_anchor_ties.
This keeps the node rows KEEP (default 0, 5, 10): every row is unchanged in place, the band's outline and plane are unchanged,
the in-between rows are deleted, and 35 x (len(KEEP)-1) quads replace the 350 (same node order pattern, so the same normals).
PM_contact is rebuilt from the new quads (same orientation). Each PM_anchor_ties row ties an anchor to a point of the PM given
by weights over old PM nodes; the point's (row, column) index position is the weighted mean of its nodes' indices (the rows
are evenly spaced along every column), and the new weights are the bilinear shape functions of the coarse quad holding it.
usage (module): t2, msg = coarsen(t, keep=(0, 5, 10))"""
import re

import numpy as np

NR, NC = 11, 36


def coarsen(t, keep=(0, 5, 10)):
    X = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(',')])
         for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t)}
    m = re.search(r'(<Elements type="quad4" name="PM_Plane">\n)(.*?)(\t*</Elements>)', t, re.S)
    E = [(int(a), [int(v) for v in b.split(',')]) for a, b in re.findall(r'<elem id="(\d+)">([^<]+)</elem>', m.group(2))]
    base = min(v for _, e in E for v in e)
    assert len({v for _, e in E for v in e}) == NR * NC and len(E) == (NR - 1) * (NC - 1)
    rc = lambda v: divmod(v - base, NC)  # noqa: E731
    # the node order pattern of an old quad, as (dr, dc) offsets from its lowest node
    e0 = E[0][1]
    r0, c0 = map(min, zip(*[rc(v) for v in e0]))
    pat = [(rc(v)[0] - r0, rc(v)[1] - c0) for v in e0]
    keep = list(keep)
    nid = lambda r, c: base + NC * r + c  # noqa: E731
    newE, eid = [], E[0][0]
    for i in range(len(keep) - 1):
        for c in range(NC - 1):
            newE.append((eid, [nid(keep[i] if dr == 0 else keep[i + 1], c + dc) for dr, dc in pat]))
            eid += 1
    ebl = ''.join(f'\t\t\t<elem id="{a}">{",".join(map(str, e))}</elem>\n' for a, e in newE)
    t = t[:m.start(2)] + ebl + t[m.end(2):]
    # element ids after the PM block must not clash: the PM block keeps its first ids, fewer elements is fine
    gone = {nid(r, c) for r in range(NR) if r not in keep for c in range(NC)}
    # PM_contact (if present): the same orientation as the old faces (old face = element node order, or reversed)
    sm = re.search(r'(<Surface name="PM_contact">\n)(.*?)(\t*</Surface>)', t, re.S)
    faces = []
    if sm:
        old_faces = [[int(v) for v in f.split(',')] for f in re.findall(r'>([\d,]+)</quad4>', sm.group(2))]
        same = any(f == E[0][1] for f in old_faces)
        rev = any(f == E[0][1][::-1] or f == [E[0][1][0]] + E[0][1][1:][::-1] for f in old_faces)
        assert same or rev, 'PM_contact orientation not recognised'
        faces = [e if same else [e[0]] + e[1:][::-1] for _, e in newE]
        t = t[:sm.start(2)] + ''.join(f'\t\t\t<quad4 id="{i + 1}">{",".join(map(str, f))}</quad4>\n'
                                        for i, f in enumerate(faces)) + t[sm.end(2):]
    # PM_anchor_ties: re-express each tied point on the coarse grid
    cm = re.search(r'<constraint name="PM_anchor_ties".*?</constraint>', t, re.S)
    blk = cm.group(0)
    moved, maxdev = 0, 0.0

    def redo(mm):
        nonlocal moved, maxdev
        rows = re.findall(r'<node id="(\d+)" bc="(\w)">([^<]+)</node>', mm.group(0))
        anc, bc = int(rows[0][0]), rows[0][1]
        w = [(int(n), -float(c)) for n, _, c in rows[1:]]
        if not any(n in gone for n, _ in w):
            return mm.group(0)
        s = sum(c for _, c in w)
        pos = [(rc(n), c / s) for n, c in w]
        R = sum(r * c for (r, _), c in pos); C = sum(cc * c for (_, cc), c in pos)
        P_old = sum(X[n] * c for n, c in w) / s
        i = max(k for k in range(len(keep) - 1) if keep[k] <= R + 1e-9) if R < keep[-1] else len(keep) - 2
        ra, rb = keep[i], keep[i + 1]
        v = (R - ra) / (rb - ra)
        c0 = min(int(np.floor(C)), NC - 2); u = C - c0
        nw = {nid(ra, c0): (1 - u) * (1 - v), nid(ra, c0 + 1): u * (1 - v), nid(rb, c0): (1 - u) * v, nid(rb, c0 + 1): u * v}
        P_new = sum(X[n] * c for n, c in nw.items())
        maxdev = max(maxdev, float(np.linalg.norm(P_new - P_old)))
        moved += 1
        ind = '\t\t\t\t'
        return ('<linear_constraint>\n' + f'{ind}<node id="{anc}" bc="{bc}">1</node>\n' +
                ''.join(f'{ind}<node id="{n}" bc="{bc}">{-c:.6f}</node>\n' for n, c in nw.items() if abs(c) > 1e-9) +
                '\t\t\t</linear_constraint>')
    nb = re.sub(r'<linear_constraint>.*?</linear_constraint>', redo, blk, flags=re.S)
    t = t.replace(blk, nb)
    # delete the dropped rows' nodes (nothing else references them: checked below)
    t = re.sub(r'\t*<node id="(\d+)">[^<]+</node>\n', lambda mm: '' if int(mm.group(1)) in gone else mm.group(0), t)
    # nothing may still point at a deleted node: elements, springs, surfaces, constraint rows, node sets
    refs = {int(v) for v in re.findall(r'\b(\d+)\b', ' '.join(re.findall(r'<(?:delem|elem|quad4|tri3)[^>]*>([^<]+)<', t)))}
    refs |= {int(v) for v in re.findall(r'<node id="(\d+)" bc=', t)}
    for ns in re.finditer(r'<NodeSet name="[^"]+">(.*?)</NodeSet>', t, re.S):
        refs |= {int(v) for v in (re.findall(r'id="(\d+)"', ns.group(1)) or re.split(r'[,\s]+', ns.group(1).strip())) if v.isdigit()}
    bad = sorted(refs & gone)
    assert not bad, f'dropped PM nodes still referenced: {bad[:10]}'
    contact = f'PM_contact rebuilt ({len(faces)} faces)' if faces else 'no PM_contact surface in this model'
    msg = (f'NOT IN SOURCE (PM shell remeshed across its width; tools/pm_coarsen.py): rows {keep} of the 11 x 36 PM grid kept '
           f'(outline, plane and the clamped row unchanged), {len(newE)} quads replace 350 (width per row from 0.36-0.68 to '
           f'{0.36 * (keep[1] - keep[0]):.1f}-{0.68 * (keep[1] - keep[0]):.1f} mm under the 2 mm shell), {len(gone)} nodes removed; '
           f'{contact}; {moved} PM_anchor_ties rows re-expressed on the coarse quads (tied point moved <= {maxdev:.3f} mm)')
    return t, msg


if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from variants7 import JOBS
    from build_batch162 import write
    RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
    src, name = sys.argv[1], sys.argv[2]
    keep = tuple(int(v) for v in sys.argv[3].split(',')) if len(sys.argv) > 3 else (0, 5, 10)
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    t2, msg = coarsen(t, keep)
    print(msg)
    write(name, t2, f'base: claude_diag/runs/{src}/{src}.feb\n{msg}\n')
