"""Batch 307 (2026-10-17; the user: "replace the AVW-Para chains + loft with solid ropes (fitted to the connector curves) on the
chains model, add rope <-> PM arch contact, and run MAX"; FEBio's spring contact cannot be set up from a .feb file, skill gotcha
33, so each connector becomes a thin solid that ordinary contact can see).
NOT IN SOURCE: every AVW-Para connector (35 a side, anchor RP 10578 / 10579 -> a tissue node) becomes a solid rope: a column of
NSEG hex8 with a S x S mm square section along the straight anchor-tissue line, its own material (solid mixture: neo-Hookean
E 0.1, nu 0 + fiber-exp-pow along the rope, fibre "local 1,2"), fitted to the connector's force-elongation table x the case
scale for the rope's rest length (ropes/rope_fit.py; F = 2 A0 lam Psi'; checked in FEBio on a single 40 mm rope: within 1-4 %
from 10 to 34 mm). The anchor-end face nodes are fixed (the anchors are fixed); the tissue-end face nodes are tied to the tissue
node (linear constraints, junction_ties settings). Removed: the AVW-Para chain springs (their <discrete> lines and discrete
sets; the discrete materials stay so dmat positions are unchanged), the chain nodes, the AVW-Para loft sheets (elements and
domains; their materials stay unused).
Contact AVW_Para_rope_PM: sliding-elastic, the ropes' side faces (primary) vs the PM arch (PM_Plane facets reversed, facing the
ropes on its -normal side as the sheets were); penalty PEN auto, two_pass 0, search_radius SR, offset 0, seg_up 2 (as walls_LA).
usage: py -3.10 build_batch307.py SRC_RUN NEW_RUN [--scale 0.5] [--s 0.5] [--nseg 10] [--pen 10] [--sr 3] [--nocontact] [--matE 0.1]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ropes'))
from junction_view_data import load  # noqa: E402
from build_batch278 import lc  # noqa: E402
from build_batch162 import write  # noqa: E402
from rope_fit import fit  # noqa: E402

SIDES = {'L': 10578, 'R': 10579}


def chains(t, side, anchor):
    ds = re.search(rf'<DiscreteSet name="AVW-Para-{side}_conn_chain">(.*?)</DiscreteSet>', t, re.S).group(1)
    pairs = [tuple(map(int, p)) for p in re.findall(r'<delem>(\d+),(\d+)</delem>', ds)]
    blk = re.search(rf'<Nodes name="AVW-Para-{side}_chain_nodes">(.*?)</Nodes>', t, re.S).group(1)
    inner = {int(v) for v in re.findall(r'<node id="(\d+)"', blk)}
    nb = {}
    for a, b in pairs:
        nb.setdefault(a, []).append(b); nb.setdefault(b, []).append(a)
    ends = []
    for first in nb[anchor]:
        prev, cur = anchor, first
        while cur in inner:
            nxt = [v for v in nb[cur] if v != prev]
            assert len(nxt) == 1, (cur, nb[cur])
            prev, cur = cur, nxt[0]
        ends.append(cur)
    assert len(ends) == 35, len(ends)
    return ends


def main(src, new, scale, s, nseg, pen, sr, contact, mate=0.1):
    t, X, E = load(src)
    nid = max(X) + 1; eid = max(E) + 1
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    rows_n, blocks, mats, doms, fixed, ties, faces, fitlog = [], [], [], [], [], '', [], []
    for side, anchor in SIDES.items():
        for k, tn in enumerate(chains(t, side, anchor), 1):
            A, T = X[anchor], X[tn]; L0 = np.linalg.norm(T - A); u = (T - A) / L0
            e1 = np.cross(u, [0, 0, 1.0]);
            if np.linalg.norm(e1) < 0.1:
                e1 = np.cross(u, [1.0, 0, 0])
            e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
            c = {}
            for i in range(nseg + 1):
                P = A + u * L0 * i / nseg
                for p in range(2):
                    for q in range(2):
                        x = P + (p - 0.5) * s * e1 + (q - 0.5) * s * e2
                        c[i, p, q] = nid; rows_n.append(f'\t\t\t<node id="{nid}">{x[0]:.6f},{x[1]:.6f},{x[2]:.6f}</node>\n'); nid += 1
            hx = [[c[i, 0, 0], c[i + 1, 0, 0], c[i + 1, 1, 0], c[i, 1, 0], c[i, 0, 1], c[i + 1, 0, 1], c[i + 1, 1, 1], c[i, 1, 1]]
                  for i in range(nseg)]
            name = f'AVW-Para-{side}_rope_{k:02d}'
            blocks.append(f'\t\t<Elements type="hex8" name="{name}">\n' +
                          ''.join(f'\t\t\t<elem id="{eid + j}">{",".join(map(str, h))}</elem>\n' for j, h in enumerate(hx)) + '\t\t</Elements>\n')
            eid += nseg
            ksi, alpha, beta, err = fit(L0, scale=scale, A0=s * s)
            fitlog.append(err)
            mats.append(f'\t\t<material id="{mid}" name="{name}_mat" type="solid mixture">\n\t\t\t<density>1.06e-09</density>\n'
                        f'\t\t\t<solid type="neo-Hookean">\n\t\t\t\t<E>{mate:g}</E>\n\t\t\t\t<v>0</v>\n\t\t\t</solid>\n'
                        f'\t\t\t<solid type="fiber-exp-pow">\n\t\t\t\t<ksi>{ksi:.6g}</ksi>\n\t\t\t\t<alpha>{alpha:.6g}</alpha>\n'
                        f'\t\t\t\t<beta>{beta:.6g}</beta>\n\t\t\t\t<fiber type="local">1,2</fiber>\n\t\t\t</solid>\n\t\t</material>\n')
            mid += 1
            doms.append(f'\t\t<SolidDomain name="{name}" mat="{name}_mat" />\n')
            fixed += [c[0, p, q] for p in range(2) for q in range(2)]
            for p in range(2):
                for q in range(2):
                    ties += lc(c[nseg, p, q], [(tn, 1.0)])
            for i in range(nseg):
                faces += [[c[i, 0, 0], c[i + 1, 0, 0], c[i + 1, 0, 1], c[i, 0, 1]], [c[i, 1, 0], c[i, 1, 1], c[i + 1, 1, 1], c[i + 1, 1, 0]],
                          [c[i, 0, 0], c[i, 1, 0], c[i + 1, 1, 0], c[i + 1, 0, 0]], [c[i, 0, 1], c[i + 1, 0, 1], c[i + 1, 1, 1], c[i, 1, 1]]]
    # remove the AVW-Para chains, chain nodes and loft sheets
    for side in SIDES:
        for pat in (rf'\t\t<Nodes name="AVW-Para-{side}_chain_nodes">.*?</Nodes>\n',
                    rf'\t\t<Elements type="tri3" name="AVW-Para-{side}_sheet">.*?</Elements>\n',
                    rf'\t\t<ShellDomain name="AVW-Para-{side}_sheet".*?</ShellDomain>\n',
                    rf'\t\t<DiscreteSet name="AVW-Para-{side}_conn_chain">.*?</DiscreteSet>\n',
                    rf'\t\t<discrete dmat="\d+" discrete_set="AVW-Para-{side}_conn_chain" />\n'):
            t, nsub = re.subn(pat, '', t, count=1, flags=re.S)
            assert nsub == 1, pat
    # add
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="AVW-Para_rope_nodes">\n' + ''.join(rows_n) + '\t\t</Nodes>\n' + t[i:]
    i = t.rindex('</Elements>') + len('</Elements>\n')
    t = t[:i] + ''.join(blocks) + t[i:]
    i = t.index('\t\t<SurfacePair name=')
    surf = ('\t\t<NodeSet name="AVW-Para_rope_anchor_ends">' + ','.join(map(str, fixed)) + '</NodeSet>\n')
    pm = [c for nm, c in E.values() if nm == 'PM_Plane']
    if contact:
        surf += ('\t\t<Surface name="AVW_Para_rope_PM_primary">\n' +
                 ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, f))}</quad4>\n' for k, f in enumerate(faces, 1)) +
                 '\t\t</Surface>\n\t\t<Surface name="AVW_Para_rope_PM_secondary">\n' +
                 ''.join(f'\t\t\t<quad4 id="{k}">{",".join(map(str, q[::-1]))}</quad4>\n' for k, q in enumerate(pm, 1)) +
                 '\t\t</Surface>\n\t\t<SurfacePair name="AVW_Para_rope_PM">\n\t\t\t<primary>AVW_Para_rope_PM_primary</primary>\n'
                 '\t\t\t<secondary>AVW_Para_rope_PM_secondary</secondary>\n\t\t</SurfacePair>\n')
    t = t[:i] + surf + t[i:]
    i = t.rindex('</material>') + len('</material>\n')
    t = t[:i] + ''.join(mats) + t[i:]
    i = t.rindex('</SolidDomain>') + len('</SolidDomain>\n') if '</SolidDomain>' in t else t.rindex('<SolidDomain ')
    i = t.index('\n', t.rindex('<SolidDomain ')) + 1
    t = t[:i] + ''.join(doms) + t[i:]
    bc = re.search(r'\t\t<bc name="BC-RP-ParaSupport" .*?</bc>\n', t, re.S).group(0)
    t = t.replace(bc, bc + bc.replace('BC-RP-ParaSupport', 'AVW-Para_rope_anchor_ends'), 1)
    jt = re.search(r'\t\t<constraint name="junction_ties" type="linear constraint">\n(.*?)\t\t\t<linear_constraint>', t, re.S)
    head = '\t\t<constraint name="AVW_Para_rope_ties" type="linear constraint">\n' + jt.group(1)
    i = t.index('</Constraints>')
    t = t[:i] + head + ties + '\t\t</constraint>\n\t' + t[i:]
    if contact:
        a = t.index('<contact name="walls_LA"'); b = t.index('</contact>', a) + len('</contact>\n')
        cblk = (t[a:b].replace('name="walls_LA" surface_pair="walls_LA"', 'name="AVW_Para_rope_PM" surface_pair="AVW_Para_rope_PM"')
                .replace('<penalty>0.5</penalty>', f'<penalty>{pen:g}</penalty>').replace('<two_pass>1</two_pass>', '<two_pass>0</two_pass>')
                .replace('<search_radius>3.0</search_radius>', f'<search_radius>{sr:g}</search_radius>').replace('<offset>2.0</offset>', '<offset>0</offset>'))
        for x in (f'<penalty>{pen:g}</penalty>', '<two_pass>0</two_pass>', f'<search_radius>{sr:g}</search_radius>', '<offset>0</offset>'):
            assert x in cblk, x
        t = t[:b] + '\t\t' + cblk + t[b:]
    msg = (f'NOT IN SOURCE: AVW-Para chains + loft sheets replaced by 70 solid ropes ({nseg} hex8, {s:g} x {s:g} mm; solid mixture '
           f'neo-Hookean E {mate:g} nu 0 + fiber-exp-pow fitted per rope to the connector table x{scale:g}, max fit error '
           f'{100 * max(fitlog):.1f} %); anchor ends fixed, tissue ends tied (AVW_Para_rope_ties); ' +
           (f'contact AVW_Para_rope_PM (rope side faces vs PM_Plane reversed; penalty {pen:g} auto, two_pass 0, search_radius {sr:g}, '
            'offset 0, seg_up 2)' if contact else 'no contact'))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch307.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--scale', 0.5), g('--s', 0.5), int(g('--nseg', 10)), g('--pen', 10), g('--sr', 3), '--nocontact' not in a, g('--matE', 0.1))
