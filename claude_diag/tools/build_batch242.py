"""Batch 242 (2026-10-09, the user, tube line): "1) Change the material of the curved edges to cervix instead of avw when they
are added to the side of the cervix (including the top cap if that worked). 2) move the attachments (CL, US, Para, etc.) from
the side of the avw/pvw/cervix to the curved sides. Attach them not at the node where the avw meets the pvw, but the node above
that." The user's choices: (1) only the U's AVW half where it lies on the cervix + the apex cap; (2) each side's own neighbour
of the seam node; test each alone and both.
  mat     canal_wrap_hex_avw hexes of the curve rings whose AVW-side leg is on the cervix (_PickedSet346) move to a new domain
          canal_wrap_hex_avw_cervix (Vagina_Cervix, the curves' damping); the apex cap apex_contour -> Vagina_Cervix.
  attach  the wall end of every CL / AVW-Para spring moves to the curve's outer surface (row 4) at column 2 of its own ring
          (column 3 = where the U's AVW half meets its PVW half; column 0 = the AVW / cervix-side leg, where they attach now);
          USL / Parcus springs move to column 4 (column 6 = the PVW-side leg). Springs on the cervix beyond the curves' end go
          to the apex cap's ring nearest in depth. Spring ends not on the curves or the cap (Parcus at the PeB) are kept.
          The spring tables are unchanged (force vs. change of length from the new rest length).
usage: py -3.10 build_batch242.py SRC NAME mat|attach|mat,attach"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from curve_rings import rings  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
COL = {'CL-': 2, 'AVW-Para': 2, 'USL-': 4, 'Parcus': 4}
LABEL = 'NOT IN SOURCE (the user, 2026-10-09, tube line)'


def cap_rings(t, grids):
    """the apex cap's rings, walked up from the curves' end cap through apex_contour's hexes"""
    mesh = ET.fromstring(t).find('Mesh')
    H = [[int(v) for v in el.text.split(',')] for e in mesh.findall('Elements') if e.get('name') == 'apex_contour' for el in e]
    byface = defaultdict(list)
    for h in H:
        byface[frozenset(h[:4])].append(h); byface[frozenset(h[4:])].append(h)
    out, used = [grids[0]], set()
    while True:
        G = out[-1]; N = -np.ones_like(G)
        for a in range(4):
            for b in range(6):
                q = [G[a, b], G[a, b + 1], G[a + 1, b + 1], G[a + 1, b]]
                hs = [h for h in byface[frozenset(q)] if tuple(h) not in used]
                if not hs:
                    continue
                h = hs[0]; used.add(tuple(h))
                for v in q:
                    N[np.where(G == v)] = h[(h.index(v) + 4) % 8]
        if (N < 0).any():
            return out[1:]
        out.append(N)


def main(src, name, mods):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    dom = defaultdict(set)
    for e in mesh.findall('Elements'):
        for el in e:
            for v in el.text.split(','):
                dom[int(v)].add(e.get('name'))
    R = rings(t)
    log = []
    if 'mat' in mods:
        ring_on_cx = {}
        for sgn, (grids, c0, B, X) in R.items():
            for k, g in enumerate(grids):
                ring_on_cx[(sgn, k)] = any('_PickedSet346' in dom[int(v)] for v in g[:, 0])
        cxn = set()
        for sgn, (grids, c0, B, X) in R.items():
            ks = [k for k in range(len(grids)) if ring_on_cx[(sgn, k)]]
            assert ks == list(range(len(ks))), ks          # one run from the apical end
            for k in ks:
                cxn |= {int(v) for v in grids[k].ravel()}
        blk = re.search(r'(\t\t<Elements type="hex8" name="canal_wrap_hex_avw">\n)(.*?)(\t\t</Elements>\n)', t, re.S)
        keep, move = [], []
        for line in blk.group(2).splitlines(keepends=True):
            ids = [int(v) for v in re.search(r'>([\d,]+)<', line).group(1).split(',')]
            (move if set(ids) <= cxn else keep).append(line)
        t = t.replace(blk.group(0), blk.group(1) + ''.join(keep) + blk.group(3)
                      + '\t\t<Elements type="hex8" name="canal_wrap_hex_avw_cervix">\n' + ''.join(move) + blk.group(3))
        t = t.replace('\t\t<SolidDomain name="canal_wrap_hex_avw" mat="Vagina_AVW" />\n',
                      '\t\t<SolidDomain name="canal_wrap_hex_avw" mat="Vagina_AVW" />\n'
                      '\t\t<SolidDomain name="canal_wrap_hex_avw_cervix" mat="Vagina_Cervix" />\n')
        assert '<SolidDomain name="apex_contour" mat="Vagina_AVW" />' in t
        t = t.replace('<SolidDomain name="apex_contour" mat="Vagina_AVW" />', '<SolidDomain name="apex_contour" mat="Vagina_Cervix" />')
        cb = re.search(r'<body_load name="wall_damping_canal_wrap_hex_avw"[^>]*>\s*<C>([^<]+)</C>', t)
        i = t.index('\t</Loads>')
        t = t[:i] + (f'\t\t<body_load name="wall_damping_canal_wrap_hex_avw_cervix" type="mass damping" '
                     f'elem_set="canal_wrap_hex_avw_cervix">\n\t\t\t<C>{cb.group(1)}</C>\n\t\t</body_load>\n') + t[i:]
        nr = {sgn: sum(ring_on_cx[(sgn, k)] for k in range(len(R[sgn][0]))) for sgn in R}
        log.append(f'{LABEL}: curves beside the cervix -> Vagina_Cervix: {len(move)} canal_wrap_hex_avw hexes (the U\'s AVW half, '
                   f'rings on the cervix: right {nr[1]}, left {nr[-1]} from the apical end) -> canal_wrap_hex_avw_cervix '
                   f'(Vagina_Cervix, damping C {cb.group(1)}); {len(keep)} stay Vagina_AVW; apex cap apex_contour -> Vagina_Cervix')
    if 'attach' in mods:
        where, caps = {}, {}
        for sgn, (grids, c0, B, X) in R.items():
            for k, g in enumerate(grids):
                for a in range(5):
                    for b in range(7):
                        where[int(g[a, b])] = (sgn, k)
            caps[sgn] = cap_rings(t, grids)
        X = R[1][3]; c0, B = R[1][1], R[1][2]
        depth = lambda v: ((X[v] - c0) @ B)[2]  # noqa: E731
        moved, kept = Counter(), Counter()
        for m in re.finditer(r'<DiscreteSet name="([^"]+)">(.*?)</DiscreteSet>', t, re.S):
            fam = m.group(1)
            col = next((c for p, c in COL.items() if fam.startswith(p)), None)
            if col is None:
                continue
            body = m.group(2)

            def rep(mm):
                a, w = int(mm.group(1)), int(mm.group(2))
                if w in where:
                    sgn, k = where[w]
                    nw = int(R[sgn][0][k][4, col])
                else:
                    sgn = 1 if X[w][0] > 0 else -1
                    cr = caps[sgn]
                    if not ({'_PickedSet346', '_PickedSet64', '_PickedSet347'} & dom[w]) or not cr or depth(w) < depth(int(cr[0][2, 3])) - 1.5:
                        kept[fam] += 1
                        return mm.group(0)
                    k = int(np.argmin([abs(depth(int(g[2, 3])) - depth(w)) for g in cr]))
                    nw = int(cr[k][4, col])
                moved[fam] += 1
                return f'<delem>{a},{nw}</delem>'
            nb = re.sub(r'<delem>(\d+),(\d+)</delem>', rep, body)
            t = t.replace(m.group(0), m.group(0).replace(body, nb))
        log.append(f'{LABEL}: spring wall ends moved to the curves\' outer surface (row 4): CL / AVW-Para -> column 2, USL / Parcus '
                   f'-> column 4 (column 3 = the AVW / PVW halves\' meeting node); moved {dict(moved)}; kept {dict(kept)}')
    print('\n'.join(log))
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch242.py)\n' + '\n'.join(log) + '\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3].split(','))
