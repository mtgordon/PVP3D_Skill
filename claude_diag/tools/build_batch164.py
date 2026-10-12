"""Batch 164 (2026-10-05, the user: "make a version of p1 with the vaginal wall being 30% less stiff and can you add contact
between the outside of the curved edge and itself? It should start to fold back on itself like a sock").
From L163_thin275_wrapend_avwbend_full_fixed (the 4 PM_avw_bottom springs on the curves' ends, PM-LA law at full, PM ends
pinned as the source):
  * self-contact of the curves' outer surface (NOT IN SOURCE): every boundary face of canal_wrap_hex_avw / _pvw that is not in
    the lumen contact (SlidingElastic1), not an end cap (normal within ~45 deg of the local canal axis, from the apex line), not
    at a curve root (all 4 nodes tied to (canal_wrap_ends) or shared with a wall) and has no node embedded in the PVW / PeB
    (canal_wrap_embed); both sides and both halves in one surface (curve_outer), used as primary AND secondary (FEBio's
    self-contact: the same surface on both sides; developers' forum advice: two_pass off, a small search_radius). Settings as
    the lumen contact SlidingElastic1 (sliding-elastic, penalty 5, auto_penalty, seg_up 2) but two_pass 0, search_radius 0.5.
  * the walls 30 % less stiff (NOT FAITHFUL): build_batch145 walls=0.7 (Vagina_AVW, Vagina_Cervix, Vagina_PVW Yeoh c1, c2,
    k x 0.7; the curves use Vagina_AVW / Vagina_PVW), then the paper's P1.
  L164_avwbend_fixed_self                 L163 + the self-contact (healthy, no wall change; the base of the next one)
  L164_avwbend_fixed_walls07_P1           L163 + walls x0.7 + P1                     (one change from L163_..._P1)
  L164_avwbend_fixed_self_walls07_P1      L163 + self-contact + walls x0.7 + P1     (one change from the one above)
usage: py -3.10 build_batch164.py"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L163_thin275_wrapend_avwbend_full_fixed'
SELF = 'L164_avwbend_fixed_self'
CURVES = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
WALLS = ('_PickedSet347', '_PickedSet346', '_PickedSet64', '_PickedSet66')
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]   # FEBio hex8 faces, outward
SEARCH_RADIUS = 0.5


def elems(t, name):
    m = re.search(r'<Elements [^>]*name="%s"[^>]*>(.*?)</Elements>' % re.escape(name), t, re.S)
    out = [[int(v) for v in c.split(',')] for c in re.findall(r'<elem id="\d+">([^<]+)</elem>', m.group(1))]
    # 2026-10-06: a local zone domain <name>_dist (build_batch194) belongs to <name>
    d = re.search(r'<Elements [^>]*name="%s_dist"[^>]*>(.*?)</Elements>' % re.escape(name), t, re.S)
    if d:
        out += [[int(v) for v in c.split(',')] for c in re.findall(r'<elem id="\d+">([^<]+)</elem>', d.group(1))]
    return out


def surf(t, name):
    i = t.index(f'<Surface name="{name}">')
    blk = t[i:t.index('</Surface>', i)]
    return {frozenset(int(v) for v in m.split(',')) for m in re.findall(r'>([\d,]+)</(?:quad4|tri3)>', blk)}


def constrained(t, name):
    i = t.index(f'<constraint name="{name}"')
    return {int(v) for v in re.findall(r'<node id="(\d+)" bc="x">1<', t[i:t.index('</constraint>', i)])}


def outer_faces(t):
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    E = {d: elems(t, d) for d in CURVES + WALLS}
    cnt, owner = Counter(), {}
    for d, es in E.items():
        for e in es:
            for f in HF:
                k = frozenset(e[i] for i in f)
                cnt[k] += 1
                if d in CURVES:
                    owner[k] = (d, tuple(e[i] for i in f), e)
    lumen = surf(t, 'SlidingElastic1Primary') | surf(t, 'SlidingElastic1Secondary')
    emb = constrained(t, 'canal_wrap_embed')
    tied = constrained(t, 'canal_wrap_ends') | {n for d in WALLS for e in E[d] for n in e}
    apex = np.array([X[n] for n in {n for e in E[CURVES[0]] for n in e} & {n for e in E[CURVES[1]] for n in e}])
    cats, out = Counter(), []
    for k, c in cnt.items():
        if c != 1 or k not in owner:
            continue
        d, f, e = owner[k]
        P = np.array([X[n] for n in f])
        nrm = np.cross(P[2] - P[0], P[3] - P[1])
        nrm /= np.linalg.norm(nrm)
        hc = np.array([X[n] for n in e]).mean(0)
        if nrm @ (P.mean(0) - hc) < 0:                       # keep the facet wound outward from its hex
            f, nrm = f[::-1], -nrm
        near = apex[np.argsort(np.linalg.norm(apex - P.mean(0), axis=1))[:60]]
        ax = np.linalg.svd(near - near.mean(0))[2][0]
        if k in lumen:
            cats['lumen (inner)'] += 1
        elif k & emb:
            cats['with a node embedded in the PVW / PeB'] += 1
        elif abs(nrm @ ax) > 0.7:
            cats['end cap'] += 1
        elif all(n in tied for n in f):
            cats['curve root (tied to / shared with a wall)'] += 1
        else:
            cats['outer ' + d] += 1
            out.append(f)
    return out, cats


def add_self_contact(t):
    faces, cats = outer_faces(t)
    srf = ('\t\t<Surface name="curve_outer">\n'
           + ''.join(f'\t\t\t<quad4 id="{i + 1}">{",".join(map(str, f))}</quad4>\n' for i, f in enumerate(faces))
           + '\t\t</Surface>\n')
    i = t.index('\t\t<SurfacePair name="SlidingElastic1">')
    t = t[:i] + srf + t[i:]
    j = t.index('</SurfacePair>', t.index('<SurfacePair name="SlidingElastic1">')) + len('</SurfacePair>')
    t = (t[:j] + '\n\t\t<SurfacePair name="curve_outer_self">\n\t\t\t<primary>curve_outer</primary>\n'
         '\t\t\t<secondary>curve_outer</secondary>\n\t\t</SurfacePair>' + t[j:])
    a = t.index('<contact name="SlidingElastic1"')
    b = t.index('</contact>', a) + len('</contact>')
    con = t[a:b].replace('name="SlidingElastic1" surface_pair="SlidingElastic1"',
                         'name="curve_outer_self" surface_pair="curve_outer_self"')
    con, k1 = re.subn(r'<two_pass>\d</two_pass>', '<two_pass>0</two_pass>', con)
    con, k2 = re.subn(r'<search_radius>[^<]*</search_radius>', f'<search_radius>{SEARCH_RADIUS}</search_radius>', con)
    assert k1 == k2 == 1
    c = t.index('<contact name="walls_LA"')
    c = t.index('</contact>', c) + len('</contact>')
    t = t[:c] + '\n\t\t' + con + t[c:]
    msg = (f'NOT IN SOURCE: self-contact of the curves\' outer surface (contact curve_outer_self, surface curve_outer as primary '
           f'and secondary, {len(faces)} faces; the curves\' boundary faces by kind: {dict(cats)}); sliding-elastic as '
           f'SlidingElastic1 but two_pass 0, search_radius {SEARCH_RADIUS}')
    return t, msg


if __name__ == '__main__':
    p = os.path.join(RUNS, SRC, SRC + '.feb')
    t, msg = add_self_contact(open(p, encoding='utf-8').read())
    print(msg)
    write(SELF, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch164.py)\n{msg}\n')
    print(' '.join(build(164, 'avwbend_fixed_walls07', {'walls': 0.7}, ['P1'], SRC)))
    print(' '.join(build(164, 'avwbend_fixed_self_walls07', {'walls': 0.7}, ['P1'], SELF)))
