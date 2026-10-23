"""Batch 320 (2026-10-18): build_batch312 (solid lofts + contact) for ANY lofts against ANY target shells (the loft-contacts plan
after AVW-Para: loft_contacts/loft_survey.py on L294b MAX found P-arcus through the levator (4-5 nodes a side at t 1) and the
fan (PM_avw_bottom) through the PM arch (1-3)).
NOT IN SOURCE: the loft sheets SHEETS (tri3 shells) become solid lofts (one penta6 layer TH thick, bottom = the sheet nodes,
top = new nodes along the averaged node normal; twisted wedges kept as shells), as build_batch312. With --contact:
TAG_solid_<targets>, sliding-elastic between every outward face of the solid lofts (primary) and the target shells' facets
(secondary), each target facet wound to face the nearest loft node at rest; penalty PEN auto, two_pass 0, search_radius SR,
offset 0 (other settings as walls_LA).
usage: py -3.10 build_batch320.py SRC_RUN NEW_RUN --sheets A,B --target T1,T2 --tag NAME [--th 0.49] [--contact] [--pen 1] [--sr 3] [--near MM] [--fac 0.8]
--near: only target facets within MM of a loft node at rest; --fac: a target facet is left out when a loft node lies behind
it (or on it, within 0.05 mm) within SR and within FAC x its size sideways (L321: 0.8 missed facet-edge nodes; LA_ICM blew up)"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, th, contact, pen, sr, sheets, target, tag, near=1e9, fac=0.8):
    t, X, E = load(src)
    nid = max(X) + 1; eid = max(E) + 1
    rows_n, blocks, doms, faces, nshell, lnodes = [], [], [], [], {}, set()
    for sh in sheets:
        tris = [c for nm, c in E.values() if nm == sh]
        assert tris, sh
        nsum = {}
        for a, b, c in tris:
            n = np.cross(X[b] - X[a], X[c] - X[a])
            for v in (a, b, c):
                nsum[v] = nsum.get(v, 0) + n
        lnodes |= set(nsum)
        top = {}
        for v, n in nsum.items():
            n = n / np.linalg.norm(n); p = X[v] + th * n
            top[v] = nid; rows_n.append(f'\t\t\t<node id="{nid}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>\n'); nid += 1
        TP = {v: X[v] + th * (n / np.linalg.norm(n)) for v, n in nsum.items()}

        def ok(a, b, c):   # all six corner tets of the wedge positive (FEBio rejects a twisted wedge at initialisation)
            A, B, C, A2, B2, C2 = X[a], X[b], X[c], TP[a], TP[b], TP[c]
            vol = lambda p0, p1, p2, p3: np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0)   # noqa: E731
            tets = [(A, B, C, A2), (B, C, A, B2), (C, A, B, C2), (A2, C2, B2, A), (B2, A2, C2, B), (C2, B2, A2, C)]
            return min(vol(*q) for q in tets) > 1e-6 * th ** 3
        shells = {tuple(f) for f in tris if not ok(*f)}
        nshell[sh] = len(shells)
        tris = [f for f in tris if tuple(f) not in shells]
        pens = [(a, b, c, top[a], top[b], top[c]) for a, b, c in tris]
        blocks.append(f'\t\t<Elements type="penta6" name="{sh}_solid">\n' +
                      ''.join(f'\t\t\t<elem id="{eid + k}">{",".join(map(str, p))}</elem>\n' for k, p in enumerate(pens)) + '\t\t</Elements>\n')
        eid += len(pens)
        mat = re.search(rf'<ShellDomain name="{sh}" mat="([^"]+)"', t).group(1)
        doms.append(f'\t\t<SolidDomain name="{sh}_solid" mat="{mat}" />\n')
        ecnt = Counter(frozenset(e) for a, b, c in tris for e in ((a, b), (b, c), (c, a)))
        for a, b, c in tris:
            faces.append(('tri3', [a, c, b]))                     # bottom (outward = -n)
            faces.append(('tri3', [top[a], top[b], top[c]]))      # top (outward = +n)
            for u, w in ((a, b), (b, c), (c, a)):
                if ecnt[frozenset((u, w))] == 1:                  # sheet boundary edge -> side face, outward
                    q = [u, w, top[w], top[u]]
                    cen = (X[a] + X[b] + X[c]) / 3
                    nq = np.cross(X[w] - X[u], (X[u] + th * (nsum[u] / np.linalg.norm(nsum[u]))) - X[u])
                    if np.dot(nq, (X[u] + X[w]) / 2 - cen) < 0:
                        q = q[::-1]
                    faces.append(('quad4', q))
        if shells:
            m_ = re.search(rf'(\t\t<Elements type="tri3" name="{sh}">\n)(.*?)(\t\t</Elements>\n)', t, re.S)
            rows = [r for r in re.findall(r'\t\t\t<elem id="\d+">[0-9,]+</elem>\n', m_.group(2))
                    if tuple(int(v) for v in re.search(r'>([0-9,]+)<', r).group(1).split(',')) in shells]
            t = t.replace(m_.group(0), m_.group(1) + ''.join(rows) + m_.group(3), 1)
        else:
            for pat in (rf'\t\t<Elements type="tri3" name="{sh}">.*?</Elements>\n', rf'\t\t<ShellDomain name="{sh}".*?</ShellDomain>\n'):
                t, k = re.subn(pat, '', t, count=1, flags=re.S); assert k == 1, pat
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + f'\t\t<Nodes name="{tag}_solid_loft_top">\n' + ''.join(rows_n) + '\t\t</Nodes>\n' + t[i:]
    i = t.rindex('</Elements>') + len('</Elements>\n')
    t = t[:i] + ''.join(blocks) + t[i:]
    i = t.index('\n', t.rindex('<SolidDomain ')) + 1
    t = t[:i] + ''.join(doms) + t[i:]
    msg = (f'NOT IN SOURCE: {"+".join(sheets)} as solid lofts: one layer of penta6, thickness {th:g} mm, extruded from each tri3 '
           f'(bottom = the sheet nodes, top = {len(rows_n)} new nodes along the averaged node normal), same materials; triangles that '
           f'cannot form a valid wedge kept as shells: {nshell}')
    if contact:
        cn = f'{tag}_solid_{"_".join(target)}'
        ln = np.array([X[v] for v in lnodes])
        sec, flip = [], 0
        for nm, q in E.values():
            if nm not in target:
                continue
            P = [X[v] for v in q]; cen = np.mean(P, 0)
            nq = np.cross(P[2] - P[0], P[3] - P[1]) if len(q) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
            j = int(np.argmin(np.linalg.norm(ln - cen, axis=1)))
            if np.dot(nq, ln[j] - cen) < 0:
                q = q[::-1]; flip += 1
            sec.append(q)
        # (L321) drop target facets with a loft node BEHIND them at rest (within SR, projecting inside the facet): the loft
        # already passes them there, and the contact would kick those nodes at the first step
        LN = np.array([X[v] for v in lnodes] + [np.array([float(c) for c in r.split(',')]) for r in re.findall(r'<node id="\d+">([^<]+)</node>', ''.join(rows_n))])
        keep, dropped, dropped_far = [], 0, 0
        for q in sec:
            P = [X[v] for v in q]; cen = np.mean(P, 0); rad = max(np.linalg.norm(p - cen) for p in P)
            nq = np.cross(P[2] - P[0], P[3] - P[1]) if len(q) == 4 else np.cross(P[1] - P[0], P[2] - P[0]); nq /= np.linalg.norm(nq)
            d = (LN - cen) @ nq; r = np.linalg.norm(LN - cen - d[:, None] * nq, axis=1)
            if np.min(np.linalg.norm(LN - cen, axis=1)) > near:
                dropped_far += 1; continue
            if np.any((d < 0.05) & (d > -sr) & (r < fac * rad)):
                dropped += 1
            else:
                keep.append(q)
        sec = keep
        assert sec, target
        gt = lambda q: 'quad4' if len(q) == 4 else 'tri3'   # noqa: E731
        surf = (f'\t\t<Surface name="{cn}_primary">\n' +
                ''.join(f'\t\t\t<{g} id="{k}">{",".join(map(str, c))}</{g}>\n' for k, (g, c) in enumerate(faces, 1)) +
                f'\t\t</Surface>\n\t\t<Surface name="{cn}_secondary">\n' +
                ''.join(f'\t\t\t<{gt(q)} id="{k}">{",".join(map(str, q))}</{gt(q)}>\n' for k, q in enumerate(sec, 1)) +
                f'\t\t</Surface>\n\t\t<SurfacePair name="{cn}">\n\t\t\t<primary>{cn}_primary</primary>\n'
                f'\t\t\t<secondary>{cn}_secondary</secondary>\n\t\t</SurfacePair>\n')
        i = t.index('\t\t<SurfacePair name=')
        t = t[:i] + surf + t[i:]
        a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
        c = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', f'name="{cn}" surface_pair="{cn}"')
             .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
             .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
        for s_ in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>'):
            assert s_ in c, s_
        t = t[:b] + '\t\t' + c + t[b:]
        msg += (f'; contact {cn}: {len(faces)} outward faces of the solid lofts vs {"+".join(target)} ({len(sec)} facets, {flip} '
                f'rewound to face the lofts, {dropped} left out: a loft node behind them at rest (within {fac:g} x the facet size), {dropped_far} more than {near:g} mm from the loft) (sliding-elastic, penalty {pen:g} auto, two_pass 0, search_radius {sr:g}, offset 0, seg_up 2)')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch320.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    gs = lambda k: tuple(a[a.index(k) + 1].split(','))   # noqa: E731
    main(a[0], a[1], g('--th', 0.49), '--contact' in a, g('--pen', 1.0), g('--sr', 3.0), gs('--sheets'), gs('--target'), a[a.index('--tag') + 1], g('--near', 1e9), g('--fac', 0.8))
