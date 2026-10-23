"""Batch 288 (2026-10-15; the user: "I don't like how the springs ... attaching to the material only pull at one node. Is there a
way to distribute the force"; picture spring_patch_view/patches_right.png approved: "build it as drawn and run MAX").
NOT IN SOURCE (the source's connectors are node-to-node too): every discrete-spring end on a curve node at rest z < ZCUT moves
to a new anchor node at the same place, tied (linear constraints, x / y / z, on displacements = a distributing coupling) to the
weighted average of the curve's boundary-face ring around that node (tent weights = face area / 4 per corner, normalised; see
spring_patch_preview.py). Spring stiffness, position and direction unchanged; its force is shared by the 7-9 patch nodes.
usage: py -3.10 build_batch288.py SRC_RUN NEW_RUN [--zcut -15] [--corner | --edge | --tipnode | --axial] [--keeptied | --tiedw W]
--tiedw W (2026-10-17 evening): spring ends on tied nodes keep W of the pull on their node, the rest of the patch shares 1 - W
(L288c: 1/2 at a tip corner; the user found L288c's tip hexes "a bit extreme"; --keeptied = W 1 collapsed the fold, L288e)
--only N1,N2,... (2026-10-19, L331): patch only these spring-end nodes (every spring on them moves to the anchor); the rest
stay on their single node"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from feb_model import Feb  # noqa: E402
from spring_patch_preview import patches  # noqa: E402
from build_batch278 import lc  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new, zcut=-15.0, corner=False, keeptied=False, tiedw=None, only=None):
    path = os.path.join(RUNS, src, src + '.feb')
    f = Feb(path); X = f.nodes
    t = open(path, encoding='utf-8').read()
    P, ends, _ = patches(f, zcut, corner=corner)
    if only:
        miss = sorted(set(only) - set(P)); assert not miss, ('--only nodes without a patch', miss)
        P = {n: v for n, v in P.items() if n in only}
        print('patched spring-end nodes (--only):', sorted(P))
    if keeptied or tiedw is not None:   # 2026-10-16 (L288e): spring ends on nodes tied onto the PeB / wall keep their single node, so their pull
        # goes into the tie as in the source model (in L288c a quarter of E29217's pull, on the tied 27588, dragged the free
        # neighbour 27587 7-11 mm out of its tip hex)
        tied = set()
        for name in ('junction_ties', 'canal_wrap_ends'):
            b = re.search(r'<constraint name="%s".*?</constraint>' % name, t, re.S)
            for lcs in re.findall(r'<linear_constraint>(.*?)</linear_constraint>', b.group(0), re.S):
                nn = re.findall(r'<node id="(\d+)" bc="x">', lcs)
                if nn:
                    tied.add(int(nn[0]))
        dropped = sorted(n for n in P if n in tied)
        if tiedw is None:
            P = {n: p for n, p in P.items() if n not in tied}
            print('spring ends on tied nodes kept on their node:', dropped)
        else:
            for n in dropped:
                rest = sum(w for v, w in P[n] if v != n)
                P[n] = [(n, tiedw)] + [(v, (1 - tiedw) * w / rest) for v, w in P[n] if v != n] if rest > 0 else [(n, 1.0)]
            print(f'spring ends on tied nodes reweighted to {tiedw:g} on their node:', dropped)
    nid = max(X) + 1
    anchor = {}
    for n in sorted(P):
        anchor[n] = nid; nid += 1
    # springs: replace the end node by its anchor, only inside DiscreteSet blocks
    nrep = 0

    def ds(m):
        nonlocal nrep
        def de(mm):
            nonlocal nrep
            a, b = (int(v) for v in mm.group(2).split(','))
            a2, b2 = anchor.get(a, a), anchor.get(b, b)
            nrep += (a2 != a) + (b2 != b)
            return f'{mm.group(1)}{a2},{b2}{mm.group(3)}'
        return re.sub(r'(<delem[^>]*>)\s*(\d+\s*,\s*\d+)\s*(</delem>)', de, m.group(0))
    t = re.sub(r'<DiscreteSet name="[^"]*">.*?</DiscreteSet>', ds, t, flags=re.S)
    assert nrep == sum(len(v) for n, v in ends.items() if n in P), (nrep, sum(len(v) for n, v in ends.items() if n in P))
    # anchor nodes after the last Nodes block
    b = t.rindex('</Nodes>\n') + len('</Nodes>\n')
    t = t[:b] + '\t\t<Nodes name="spring_patch_anchors">\n' + ''.join(
        f'\t\t\t<node id="{anchor[n]}">' + ','.join(f'{c:.9g}' for c in X[n]) + '</node>\n' for n in sorted(P)) + '\t\t</Nodes>\n' + t[b:]
    # constraints: a new block with junction_ties' settings
    blk = re.search(r'\t\t<constraint name="junction_ties".*?</constraint>\n', t, re.S).group(0)
    head = re.match(r'(\t\t<constraint name="junction_ties".*?\n)(.*?)\t\t\t<linear_constraint>', blk, re.S)
    rows = ''.join(lc(anchor[n], P[n]) for n in sorted(P))
    old = re.search(r'\t\t<constraint name="spring_patch_ties".*?\t\t</constraint>\n', t, re.S)
    if old:   # 10-22 (L423): append to the existing block; a second small block's near-zero multipliers never met the
        #       relative augmentation tolerance (~10 augmentations a step, 3x slower)
        t = t[:old.end() - len('\t\t</constraint>\n')] + rows + t[old.end() - len('\t\t</constraint>\n'):]
    else:
        t = t.replace(blk, blk + head.group(1).replace('junction_ties', 'spring_patch_ties') + head.group(2) + rows + '\t\t</constraint>\n', 1)
    fams = {}
    for n, v in ends.items():
        if n not in P:
            continue
        for nm, _ in v:
            fams[nm] = fams.get(nm, 0) + 1
    msg = (f'NOT IN SOURCE: {len(P)} discrete-spring end nodes on the curves at rest z < {zcut}' + (f' (only {sorted(P)}, --only)' if only else '') + f' ({nrep} spring ends: {fams}) moved '
           f'to anchor nodes (block spring_patch_anchors) tied to the tent-weighted average of their one-ring curve-face patch '
           f'({min(len(p) for p in P.values())}-{max(len(p) for p in P.values())} nodes; constraint spring_patch_ties, '
           f'settings as junction_ties); springs otherwise unchanged' + ('; spring ends on tied nodes (junction_ties / canal_wrap_ends) KEEP their single node' if keeptied and tiedw is None else '') + (f'; spring ends on tied nodes (junction_ties / canal_wrap_ends): {tiedw:g} of the pull on their node, {1 - tiedw:g} shared by the rest of their patch' if tiedw is not None else '') + ('; TIP CORNER ENDS SPREAD UP THE CURVE (node 1/2, next node up the canal 1/3, the one after 1/6)' if corner == 'axial' else '; TIP CORNER ENDS LEFT ON THEIR SINGLE NODE (no patch)' if corner == 'node' else '; TIP CORNERS ALONG THE EDGE: at the curve tips (ring faces meeting at > 60 deg) the node 1/2 + its 2 neighbours along the corner edge 1/4 each' if corner == 'edge' else '; CORNER-AWARE: at the curve tips (ring faces meeting at > 60 deg) only the faces within 45 deg of the spring direction' if corner else ''))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch288.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--zcut') + 1]) if '--zcut' in a else -15.0, 'axial' if '--axial' in a else 'edge' if '--edge' in a else 'node' if '--tipnode' in a else '--corner' in a, '--keeptied' in a, float(a[a.index('--tiedw') + 1]) if '--tiedw' in a else None, [int(v) for v in a[a.index('--only') + 1].split(',')] if '--only' in a else None)
