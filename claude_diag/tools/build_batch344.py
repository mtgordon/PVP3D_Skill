"""Batch 344 (2026-10-19 night; the user's go, item 1 of NEXT_SESSION_PROMPT_2026-10-19b): the AVW-Para lofts' curve ends tied
DIRECTLY onto the curve. In L342 (AVW-Para as solid fibre lofts on the springs line) the loft's curve edge is mixed: 21 nodes a
side are curve nodes themselves, 14 a side are spring-patch ANCHORS (build_batch288: a separate node at the same place, tied by
penalty linear constraints, spring_patch_ties, to the tent-weighted average of the curve-face ring around its curve node). Once
the AVW-Para springs were taken out (L342b) those anchors carry only the loft, whose stiff in-plane edge then pulls each anchor
against neighbours that are themselves curve nodes in the same ring: the penalty ties fight the loft (L342 crawls, ~48
iterations a step; L342a chains + sheet just as slow). Here every anchor used by the DOMS elements is replaced, in those
elements only, by the curve node it sits on (distance 0 at rest); anchors then referenced by nothing (no spring, element, set
or surface) are removed with their spring_patch_ties rows. Anchors still carrying other springs are left as they are.
NOT IN SOURCE (as build_batch288 / 337).
usage: py -3.10 build_batch344.py SRC_RUN NEW_RUN [--doms AVW-Para-L_sheet,AVW-Para-R_sheet,AVW-Para-L_sheet_solid,AVW-Para-R_sheet_solid]
[--sets SET1,SET2] (spring DiscreteSets whose anchor ends move onto the curve node too, e.g. the chains of L342a)"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
CURVE = ('canal_wrap_hex_avw', 'canal_wrap_hex_avw_cervix', 'canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_fine')
DOMS = ('AVW-Para-L_sheet', 'AVW-Para-R_sheet', 'AVW-Para-L_sheet_solid', 'AVW-Para-R_sheet_solid')


def main(src, new, doms=DOMS, sets=()):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    X = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    ablk = re.search(r'\t\t<Nodes name="spring_patch_anchors">\n(.*?)\t\t</Nodes>\n', t, re.S)
    anc = [int(i) for i in re.findall(r'<node id="(\d+)">', ablk.group(1))]
    curve = set()
    for d in CURVE:
        b = re.search(rf'<Elements type="[^"]+" name="{d}">(.*?)</Elements>', t, re.S).group(1)
        for c in re.findall(r'<elem id="\d+">([^<]+)</elem>', b):
            curve |= {int(v) for v in c.split(',')}
    cid = np.array(sorted(curve)); CX = np.array([X[i] for i in cid])
    onto = {}
    for a in anc:
        d = np.linalg.norm(CX - X[a], axis=1); j = int(np.argmin(d))
        if d[j] < 1e-6:
            onto[a] = int(cid[j])
    nrep, used = 0, set()

    def el(m):
        nonlocal nrep
        def one(mm):
            nonlocal nrep
            ns = [int(v) for v in mm.group(2).split(',')]
            ns2 = [onto.get(v, v) for v in ns]
            nrep += sum(a != b for a, b in zip(ns, ns2)); used.update(v for v in ns if v in onto)
            assert len(set(ns2)) == len(ns2), ('degenerate element', mm.group(0))
            return mm.group(1) + ','.join(map(str, ns2)) + mm.group(3)
        return m.group(1) + re.sub(r'(<elem id="\d+">)([^<]+)(</elem>)', one, m.group(2)) + m.group(3)
    for d in doms:
        t, k = re.subn(rf'(<Elements type="[^"]+" name="{d}">)(.*?)(</Elements>)', el, t, flags=re.S)
        assert k == 1, d
    def dset(m):
        nonlocal nrep
        def de(mm):
            nonlocal nrep
            a, b = (int(v) for v in mm.group(2).split(','))
            used.update(v for v in (a, b) if v in onto)
            a2, b2 = onto.get(a, a), onto.get(b, b)
            nrep += (a2 != a) + (b2 != b)
            return f'{mm.group(1)}{a2},{b2}{mm.group(3)}'
        return re.sub(r'(<delem[^>]*>)\s*(\d+\s*,\s*\d+)\s*(</delem>)', de, m.group(0))
    for s in sets:
        t, k = re.subn(rf'<DiscreteSet name="{s}">.*?</DiscreteSet>', dset, t, flags=re.S)
        assert k == 1, s
    # anchors still referenced anywhere (springs, other elements, sets, surfaces)?
    refs = set()
    for b in re.findall(r'<delem[^>]*>([^<]+)</delem>', t):
        refs |= {int(v) for v in b.split(',')}
    for b in re.findall(r'<elem id="\d+">([^<]+)</elem>', t):
        refs |= {int(v) for v in b.split(',')}
    for b in re.findall(r'<NodeSet name="[^"]*">([^<]*)</NodeSet>', t):
        refs |= {int(v) for v in b.replace('\n', ',').split(',') if v.strip()}
    for b in re.findall(r'<n id="(\d+)"', t):
        refs.add(int(b))
    drop = sorted(a for a in used if a not in refs)
    keep = sorted(a for a in used if a in refs)
    ds = set(drop)
    nb = ''.join(r for r in re.findall(r'\t\t\t<node id="\d+">[^<]+</node>\n', ablk.group(1)) if int(re.search(r'id="(\d+)"', r).group(1)) not in ds)
    t = t.replace(ablk.group(0), (f'\t\t<Nodes name="spring_patch_anchors">\n{nb}\t\t</Nodes>\n' if nb else ''), 1)
    cblk = re.search(r'\t\t<constraint name="spring_patch_ties".*?</constraint>\n', t, re.S).group(0)
    rows = re.findall(r'\t\t\t<linear_constraint>\n.*?\t\t\t</linear_constraint>\n', cblk, re.S)
    gone = [r for r in rows if int(re.search(r'<node id="(\d+)" bc="[xyz]">1</node>', r).group(1)) in ds]
    c2 = cblk
    for r in gone:
        c2 = c2.replace(r, '', 1)
    assert len(gone) == 3 * len(drop), (len(gone), len(drop))
    t = t.replace(cblk, c2, 1)
    msg = (f'NOT IN SOURCE: the AVW-Para lofts\' curve ends tied DIRECTLY onto the curve (tools/build_batch344.py): in {", ".join(doms)} '
           f'{len(used)} spring-patch anchors ({nrep} element references) replaced by the curve node each sits on (0 mm at rest); '
           f'{len(drop)} anchors left unused removed with their {len(gone)} spring_patch_ties rows; '
           f'{len(keep)} kept (still carry other springs: {keep}); map {dict(sorted((a, onto[a]) for a in used))}')
    print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch344.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], tuple(a[a.index('--doms') + 1].split(',')) if '--doms' in a else DOMS,
         tuple(a[a.index('--sets') + 1].split(',')) if '--sets' in a else ())
