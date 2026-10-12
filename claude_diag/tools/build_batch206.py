"""Batch 206 (2026-10-07, the user: if L203 opens the gap, "duplicate the 200 shared curve-PVW nodes and make the whole join
elastic, then build the tube default's healthy and P2"). One change on an L203 file (the curves' ties to the PVW columns
already elastic): the 200 nodes the PVW half of the curves (canal_wrap_hex_pvw, 8 of them also in canal_wrap_hex_pvw_apex)
SHARES with the PVW (_PickedSet64) are duplicated. The curve elements, and the curve-side faces of every surface, move to
the duplicates (new ids after the highest node, in their own <Nodes> block); the PVW, its surfaces, the support springs
(USL / Parcus DiscreteSets) and the constraints that use these nodes as wall-side masters (canal_wrap_ends, apex_close_ties)
keep the originals. Each duplicate follows its original through 3 rows (x, y, z) of a new linear constraint
canal_wrap_pvw_shared_soft, maxaug 0 (an elastic zero-length spring) at the same penalty as the L203 file's
canal_wrap_pvw_soft. NOT IN SOURCE.
A face is curve-side when its nodes lie in one curve element and in no PVW element; a face in both (the curve-PVW interface)
is curve-side only in the curve surfaces (CURVE_SURFS).
usage: py -3.10 build_batch206.py SRC_L203_RUN NEW_NAME"""
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
PVW = '_PickedSet64'
CURVES = ('canal_wrap_hex_pvw', 'canal_wrap_hex_pvw_apex')
CURVE_SURFS = ('Load-curves_surf', 'curve_all')


def elems(mesh, name):
    return [[int(v) for v in el.text.split(',')] for e in mesh.findall('Elements') if e.get('name') == name for el in e]


if __name__ == '__main__':
    src, name = sys.argv[1], sys.argv[2]
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    root = ET.fromstring(t)
    mesh = root.find('Mesh')
    pen = re.search(r'name="canal_wrap_pvw_soft".*?<penalty>([^<]+)</penalty>', t, re.S).group(1)
    pvw_el = [set(e) for e in elems(mesh, PVW)]
    cur_el = [set(e) for c in CURVES for e in elems(mesh, c)]
    pvw_n = set().union(*pvw_el)
    S = sorted(set().union(*cur_el) & pvw_n)
    xyz = {int(n.get('id')): n.text for n in mesh.iter('node')}
    top = max(xyz)
    dup = {s: top + 1 + i for i, s in enumerate(S)}
    rmap = lambda ids: [dup.get(i, i) for i in ids]  # noqa: E731

    # 1. the curve elements -> duplicates
    for c in CURVES:
        a = t.index(f'<Elements type="hex8" name="{c}"')
        b = t.index('</Elements>', a)
        blk = re.sub(r'(<elem id="\d+">)([^<]+)(</elem>)',
                     lambda m: m.group(1) + ','.join(map(str, rmap(int(v) for v in m.group(2).split(',')))) + m.group(3),
                     t[a:b])
        t = t[:a] + blk + t[b:]

    # 2. curve-side faces of every surface
    def in_any(fs, els):
        return any(fs <= e for e in els)
    report = []
    for s in mesh.findall('Surface'):
        nm = s.get('name')
        a = t.index(f'<Surface name="{nm}"')
        b = t.index('</Surface>', a)
        cnt = [0]

        def face(m):
            ids = [int(v) for v in m.group(2).split(',')]
            fs = set(ids)
            if not fs & set(dup):
                return m.group(0)
            curve, wall = in_any(fs, cur_el), in_any(fs, pvw_el)
            if curve and (not wall or nm in CURVE_SURFS):
                cnt[0] += 1
                return m.group(1) + ','.join(map(str, rmap(ids))) + m.group(3)
            return m.group(0)
        blk = re.sub(r'(<(?:quad4|tri3) id="\d+">)([^<]+)(</(?:quad4|tri3)>)', face, t[a:b])
        t = t[:a] + blk + t[b:]
        if cnt[0]:
            report.append(f'{nm} {cnt[0]}')

    # 3. the duplicate nodes, after the last <Nodes> block (ids must rise through the whole Mesh)
    i = t.rindex('</Nodes>') + len('</Nodes>\n')
    t = t[:i] + '\t\t<Nodes name="curve_pvw_dup_nodes">\n' + ''.join(
        f'\t\t\t<node id="{dup[s]}">{xyz[s]}</node>\n' for s in S) + '\t\t</Nodes>\n' + t[i:]

    # 4. each duplicate follows its original elastically
    rows = ''.join(f'\t\t\t<linear_constraint>\n\t\t\t\t<node id="{dup[s]}" bc="{d}">1</node>\n'
                   f'\t\t\t\t<node id="{s}" bc="{d}">-1</node>\n\t\t\t</linear_constraint>\n' for s in S for d in 'xyz')
    i = t.index('</constraint>', t.index('<constraint name="canal_wrap_pvw_soft"')) + len('</constraint>\n')
    t = t[:i] + (f'\t\t<constraint name="canal_wrap_pvw_shared_soft" type="linear constraint">\n\t\t\t<tol>0.01</tol>\n'
                 f'\t\t\t<penalty>{pen}</penalty>\n\t\t\t<maxaug>0</maxaug>\n' + rows + '\t\t</constraint>\n') + t[i:]

    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch206.py)\n'
                   f'NOT IN SOURCE: the {len(S)} nodes the PVW half of the curves shared with the PVW duplicated '
                   f'(ids {top + 1}-{top + len(S)}); curve elements + curve-side faces moved to them ({"; ".join(report)}); '
                   f'each duplicate tied to its original by canal_wrap_pvw_shared_soft, penalty {pen} N/mm per dof, maxaug 0; '
                   f'PVW, its surfaces, the USL / Parcus springs and the wall-side masters keep the originals\n')
    print(len(S), 'duplicated;', '; '.join(report))
