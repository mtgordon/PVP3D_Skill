"""Batch 169 (2026-10-05, the user: on L166_wrapfull_endspr_walls035_P1 "the contact between the curved walls fails around 0.78
and it penetrates itself" -> "build ... the contact fix").
Found (scratchpad self_pen.py): from t 0.78 the curves' distal ends fold back along the canal (the "sock") and 6 -> 26 of
their outer-surface nodes pass into the curve 5-10 mm further up, through its lumen-side faces and end caps, which were not in
the self-contact (curve_outer held only the outer faces). FEBio's sliding-elastic search_radius is a length (FEBio source:
ADD_PARAMETER(m_srad, "search_radius")->setUnits(UNIT_LENGTH); a node whose gap exceeds it is released), and 0.5 mm is
less than a fold step against the curves' 1.6-2.8 mm facets.
From L166_wrapfull_endspr_self_damp (the curves to the end, the end springs, self-contact, wall damping):
  * the self-contact surface (renamed curve_all) = every boundary face of the curves except their roots (all 4 nodes
    tied to / shared with a wall) and faces with a node embedded in the PVW / PeB: the outer surface, the end caps and the
    lumen-side faces (the lumen faces stay in SlidingElastic1 as well), wound outward from their hex;
  * search_radius kept at 0.5 mm: read tests of the whole-boundary surface (scratchpad readtest4/mk_tests.py) failed or crawled
    in the first step at 0.75 / 1.0 / 1.5 mm (at rest, nodes take the facets meeting them at the end-cap corners and across
    the tight inside of the U), and ran normally at 0.5; the outer faces alone ran at 1.5. With every boundary face in the
    surface, a folding node meets a facet when it first touches, so it need not go 0.5 mm in before it is caught.
Then walls x0.35 and P1 / P2 (build_batch145), as L166.
  L169_wrapfull_endspr_selfall_damp                 (healthy, source wall stiffness)
  L169_wrapfull_endspr_selfall_walls035_{P1,P2}
usage: py -3.10 build_batch169.py"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write, block  # noqa: E402
from build_batch164 import elems, constrained, HF, CURVES, WALLS  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L166_wrapfull_endspr_self_damp'
NAME = 'L169_wrapfull_endspr_selfall_damp'
SEARCH_RADIUS = 0.5


def all_faces(t):
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    E = {d: elems(t, d) for d in CURVES + WALLS}
    cnt, owner = Counter(), {}
    for d, es in E.items():
        for e in es:
            for f in HF:
                k = frozenset(e[i] for i in f)
                cnt[k] += 1
                if d in CURVES:
                    owner[k] = (tuple(e[i] for i in f), e)
    emb = constrained(t, 'canal_wrap_embed')
    tied = constrained(t, 'canal_wrap_ends') | {n for d in WALLS for e in E[d] for n in e}
    out, cats = [], Counter()
    for k, c in cnt.items():
        if c != 1 or k not in owner:
            continue
        f, e = owner[k]
        if k & emb:
            cats['embedded (left out)'] += 1
            continue
        if all(n in tied for n in f):
            cats['root (left out)'] += 1
            continue
        P = np.array([X[n] for n in f])
        nrm = np.cross(P[2] - P[0], P[3] - P[1])
        if nrm @ (P.mean(0) - np.array([X[n] for n in e]).mean(0)) < 0:
            f = f[::-1]
        out.append(f)
        cats['in'] += 1
    return out, cats


if __name__ == '__main__':
    p = os.path.join(RUNS, SRC, SRC + '.feb')
    t = open(p, encoding='utf-8').read()
    faces, cats = all_faces(t)
    i, j = block(t, 'Surface', 'curve_outer')
    t = (t[:i] + '<Surface name="curve_all">\n'
         + ''.join(f'\t\t\t<quad4 id="{q + 1}">{",".join(map(str, f))}</quad4>\n' for q, f in enumerate(faces)) + '\t\t</Surface>'
         + t[j:])
    for a, b in (('<primary>curve_outer</primary>', '<primary>curve_all</primary>'),
                 ('<secondary>curve_outer</secondary>', '<secondary>curve_all</secondary>'),
                 ('<SurfacePair name="curve_outer_self">', '<SurfacePair name="curve_all_self">'),
                 ('name="curve_outer_self" surface_pair="curve_outer_self"', 'name="curve_all_self" surface_pair="curve_all_self"')):
        assert t.count(a) == 1, a
        t = t.replace(a, b)
    a = t.index('<contact name="curve_all_self"')
    b = t.index('</contact>', a)
    con, k = re.subn(r'<search_radius>[^<]*</search_radius>', f'<search_radius>{SEARCH_RADIUS}</search_radius>', t[a:b])
    assert k == 1
    t = t[:a] + con + t[b:]
    msg = (f'NOT IN SOURCE: the curves\' self-contact on their whole boundary (surface curve_all, {len(faces)} faces: '
           f'{dict(cats)}; the outer surface, the end caps and the lumen-side faces), search_radius {SEARCH_RADIUS} mm (was 0.5; '
           f'the outer surface only)')
    print(msg)
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch169.py)\n{msg}\n')
    print(' '.join(build(169, 'wrapfull_endspr_selfall_walls035', {'walls': 0.35}, ['P1', 'P2'], NAME)))
