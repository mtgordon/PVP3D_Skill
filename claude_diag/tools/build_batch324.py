"""Batch 324 (2026-10-18 night). L312b (solid AVW-Para lofts <-> PM arch) holds the lofts far better than any shell version but
crawls from t 0.637: the right loft's chain nodes (29632-34) chatter on PM edge node 24199 as they slide over the PM's FREE
inner edge (the L296 crawl). L296e's flat 3 mm ledge only moved the chatter to the ledge's own edge. Here (NOT IN SOURCE, a
convergence aid like the ledges): a CURLED lip along every PM boundary edge off the clamped arc, NR rows of SEG mm quads
turning away from the lofts' side (toward +n of the PM facets, the lofts lie on -n) by ANG degrees in all, so a loft sliding
off the edge rolls round a rounded surface and never reaches a free end. PM_Yeoh, its own shell domain PM_lip, thickness TH
(thinner than the PM's 2 mm so the curl is not pinched); its facets join the loft <-> PM contact's secondary surface, wound
as the PM's (facing the lofts). Unloaded, attached to nothing else.
usage: py -3.10 build_batch324.py SRC_RUN NEW_RUN [--rows 4] [--seg 1.0] [--ang 120] [--th 0.5] [--surf AVW_Para_solid_PM_secondary]"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, nr, seg, ang, th, surf):
    t, X, E = load(src)
    pm = [c for nm, c in E.values() if nm == 'PM_Plane']
    fixed = {int(v) for v in re.findall(r'\d+', re.search(r'<NodeSet name="BC-PM-outer-arc">(.*?)</NodeSet>', t, re.S).group(1))}
    cnt = Counter(frozenset((q[i], q[(i + 1) % 4])) for q in pm for i in range(4))
    edges = []
    for q in pm:
        for i in range(4):
            a, b = q[i], q[(i + 1) % 4]
            if cnt[frozenset((a, b))] == 1 and not (a in fixed and b in fixed):
                edges.append((a, b, q))
    out, nrm = defaultdict(list), defaultdict(list)
    for a, b, q in edges:
        P = np.array([X[v] for v in q]); n = np.cross(P[2] - P[0], P[3] - P[1]); n /= np.linalg.norm(n)
        o = np.cross(X[b] - X[a], n); o /= np.linalg.norm(o)
        if np.dot(o, (X[a] + X[b]) / 2 - P.mean(0)) < 0:
            o = -o
        for v in (a, b):
            out[v].append(o); nrm[v].append(n)
    # which side are the lofts on? (expected -n)
    sec = re.search(rf'<Surface name="{surf}">(.*?)</Surface>', t, re.S)
    assert sec, surf
    nid = max(X) + 1; eid = max(E) + 1
    ring = {v: [v] for v in out}; rows_n = []
    for v in sorted(out):
        o = np.mean(out[v], 0); o /= np.linalg.norm(o); n = np.mean(nrm[v], 0); n -= np.dot(n, o) * o; n /= np.linalg.norm(n)
        p = X[v].copy()
        for k in range(1, nr + 1):
            th_k = np.radians(ang) * (k - 0.5) / nr
            p = p + seg * (np.cos(th_k) * o + np.sin(th_k) * n)
            ring[v].append(nid); rows_n.append(f'\t\t\t<node id="{nid}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>\n'); nid += 1
    quads = [(ring[b][k], ring[a][k], ring[a][k + 1], ring[b][k + 1]) for a, b, q in edges for k in range(nr)]
    rows_e = [f'\t\t\t<elem id="{eid + k}">{",".join(map(str, c))}</elem>\n' for k, c in enumerate(quads)]
    # rest clearance from every loft node (primary of the contact) to the new lip nodes
    m = re.search(rf'<Surface name="{surf.replace("_secondary", "_primary")}">(.*?)</Surface>', t, re.S)
    lp = np.array([X[int(v)] for r in re.findall(r'>([0-9,]+)<', m.group(1)) for v in r.split(',') if int(v) in X]) if m else np.zeros((0, 3))
    newp = np.array([[float(c) for c in re.search(r'>([^<]+)<', r).group(1).split(',')] for r in rows_n])
    clear = min(np.min(np.linalg.norm(lp - p, axis=1)) for p in newp) if len(lp) else float('nan')
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="PM_lip_nodes">\n' + ''.join(rows_n) + '\t\t</Nodes>\n' + t[i:]
    i = t.rindex('</Elements>') + len('</Elements>\n')
    t = t[:i] + '\t\t<Elements type="quad4" name="PM_lip">\n' + ''.join(rows_e) + '\t\t</Elements>\n' + t[i:]
    dom = re.search(r'\t\t<ShellDomain name="PM_Plane".*?</ShellDomain>\n', t, re.S).group(0)
    d2 = re.sub(r'<shell_thickness>[^<]*</shell_thickness>', f'<shell_thickness>{th:g}</shell_thickness>', dom.replace('name="PM_Plane"', 'name="PM_lip"'))
    t = t.replace(dom, dom + d2, 1)
    a = t.index(f'<Surface name="{surf}">'); b = t.index('\t\t</Surface>', a)
    tag = re.findall(r'<(quad4|tri3) id=', t[a:b])
    k0 = len(tag)
    add = ''.join(f'\t\t\t<quad4 id="{k0 + k + 1}">{",".join(map(str, c[::-1]))}</quad4>\n' for k, c in enumerate(quads))
    t = t[:b] + add + t[b:]
    msg = (f'NOT IN SOURCE (convergence aid): curled PM lip, {len(quads)} quad4 (PM_lip, PM_Yeoh, thickness {th:g}) in {nr} rows of '
           f'{seg:g} mm along the {len(edges)} PM boundary edges off the clamped arc, turning {ang:g} deg away from the lofts (+n); '
           f'added to {surf} (wound as the PM facets); nearest contact-primary node to a lip node at rest {clear:.2f} mm')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch324.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], int(g('--rows', 4)), g('--seg', 1.0), g('--ang', 120), g('--th', 0.5),
         a[a.index('--surf') + 1] if '--surf' in a else 'AVW_Para_solid_PM_secondary')
