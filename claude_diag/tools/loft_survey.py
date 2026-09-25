"""Survey every loft (ShellDomain *_fan) against the Abaqus CONN3D2 family it stands in for.

usage: py -3.10 loft_survey.py [MODEL.feb]      (default: runs/LPFmsTk_DM1/LPFmsTk_DM1.feb)
Per family: connector count and behaviours, where each end is in the FEBio mesh (by coordinate), whether the
ends are loft nodes, connector length range, loft area / node count, the length of the loft's two connector-end
edges, and what holds the loft's anchor-side nodes (BCs, other domains). Also imported by loft_fit_all.py.
"""
import os
import re
import sys
from collections import defaultdict, Counter

import numpy as np

import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import INP_FILE, RUNS_DIR
from febmodel import Feb

FAMILY_LOFT = {   # Abaqus connector family (elset name minus the trailing number) -> FEBio loft domain
    'P-arcus-L': 'P-arcus-L_fan', 'P-arcus-R': 'P-arcus-R_fan',
    'CL-L': 'CL-L_fan', 'CL-R': 'CL-R_fan', 'USL-L': 'USL-L_fan', 'USL-R': 'USL-R_fan',
    'AVW-Para-L': 'AVW-Para-L_fan', 'AVW-Para-R': 'AVW-Para-R_fan',
    'PM': 'PM_fan', 'PM-middle': 'PM_fan',
    'PM_PeB_Left_': 'PM_PeB_Left_fan', 'PM_PeB_Right_': 'PM_PeB_Right_fan',
    'PM_avw_bottom_left_': 'PM_avw_bottom_left_fan', 'PM_avw_bottom_right_': 'PM_avw_bottom_right_fan',
    'PeB-constrin': 'PeB-constrin_fan',
}


def unq(s):
    return s.strip().strip('"')


def parse_inp(path=INP_FILE):
    """Parts' nodes, instances, assembly nodes, CONN3D2 connectors (with elset + behaviour) and behaviours."""
    L = open(path, encoding='latin-1').read().splitlines()
    parts, inst, asm_nodes, conns, behav = {}, {}, {}, [], {}
    elsets = defaultdict(list)
    i, n = 0, len(L)
    part = None
    in_asm = in_inst = False
    mode = None
    cur_b = None
    while i < n:
        s = L[i].strip()
        low = s.lower()
        if low.startswith('**') or not s:
            i += 1
            continue
        if s.startswith('*'):
            mode = None
            if low.startswith('*part,'):
                part = unq(re.search(r'name=(".*?"|[^,]+)', s, re.I).group(1))
                parts[part] = {}
            elif low.startswith('*end part'):
                part = None
            elif low.startswith('*assembly'):
                in_asm = True
            elif low.startswith('*end assembly'):
                in_asm = False
            elif low.startswith('*instance,'):
                nm = unq(re.search(r'name=(".*?"|[^,]+)', s, re.I).group(1))
                pt = unq(re.search(r'part=(".*?"|[^,]+)', s, re.I).group(1))
                tr = []
                j = i + 1
                while not L[j].strip().startswith('*'):
                    tr.append([float(v) for v in L[j].split(',') if v.strip()])
                    j += 1
                inst[nm] = (pt, tr)
                in_inst = True
            elif low.startswith('*end instance'):
                in_inst = False
            elif low.startswith('*node') and not low.startswith('*node output') and not low.startswith('*node print'):
                mode = 'pnode' if part else ('anode' if in_asm and not in_inst else None)
            elif low.startswith('*element') and 'conn3d2' in low:
                mode = 'conn'
            elif low.startswith('*connector section'):
                es = unq(re.search(r'elset=(".*?"|[^,]+)', s, re.I).group(1))
                bh = unq(re.search(r'behavior=(".*?"|[^,]+)', s, re.I).group(1))
                if conns and conns[-1].get('elset') is None:
                    conns[-1]['elset'], conns[-1]['behavior'] = es, bh
            elif low.startswith('*elset,') and in_asm:
                es = unq(re.search(r'elset=(".*?"|[^,]+)', s, re.I).group(1))
                mode = ('elset', es, 'generate' in low)
            elif low.startswith('*connector behavior'):
                cur_b = unq(re.search(r'name=(".*?"|[^,]+)', s, re.I).group(1))
                behav[cur_b] = []
            elif cur_b and low.startswith('*connector'):
                behav[cur_b].append({'kw': s, 'table': []})
                mode = 'btable'
            elif low.startswith('*step'):
                cur_b = None
            i += 1
            continue
        # data lines
        v = [t.strip() for t in s.split(',') if t.strip()]
        if mode == 'pnode':
            parts[part][int(v[0])] = np.array([float(x) for x in v[1:4]])
        elif mode == 'anode':
            asm_nodes[int(v[0])] = np.array([float(x) for x in v[1:4]])
        elif mode == 'conn':
            conns.append({'eid': int(v[0]), 'a': v[1], 'b': v[2], 'elset': None, 'behavior': None})
        elif isinstance(mode, tuple) and mode[0] == 'elset':
            if mode[2]:
                a, b, st = (int(x) for x in v[:3])
                elsets[mode[1]].extend(range(a, b + 1, st))
            else:
                elsets[mode[1]].extend(int(x) for x in v if x.lstrip('-').isdigit())
        elif mode == 'btable':
            behav[cur_b][-1]['table'].append([float(x) for x in v])
        i += 1
    for c in conns:   # check the adjacency-based elset against the *Elset definitions
        if c['elset'] and elsets.get(c['elset']) and c['eid'] not in elsets[c['elset']]:
            c['elset_mismatch'] = True
    return parts, inst, asm_nodes, conns, behav


def end_xyz(ref, parts, inst, asm_nodes):
    """'Inst.node' or an assembly node number -> coordinates (instances here have no transform)."""
    if '.' in ref:
        nm, nd = ref.rsplit('.', 1)
        pt, tr = inst[unq(nm)]
        assert not tr, f'instance {nm} has a transform {tr}; not handled'
        return parts[pt][int(nd)], unq(nm)
    return asm_nodes[int(ref)], 'assembly'


def family(elset):
    return re.sub(r'-?\d+$', '', elset) if elset not in FAMILY_LOFT else elset


class NodeIndex:
    def __init__(self, feb):
        self.ids = np.array(sorted(feb.nodes))
        self.P = np.array([feb.nodes[i] for i in self.ids])

    def nearest(self, x):
        d = np.linalg.norm(self.P - x, axis=1)
        k = int(np.argmin(d))
        return int(self.ids[k]), float(d[k])


def bc_nodes(feb):
    """node id -> list of BC names that fix/prescribe it (active in the model)."""
    out = defaultdict(list)
    bnd = feb.root.find('Boundary')
    if bnd is None:
        return out
    for bc in bnd:
        ns = bc.get('node_set')
        if ns and ns in feb.nodesets:
            for nd in feb.nodesets[ns]:
                out[nd].append(bc.get('name') or bc.get('type'))
    return out


def loft_elements(feb, dom):
    return feb.elem_blocks[dom][1]


def tri_area(feb, conn):
    P = [feb.nodes[n] for n in conn]
    if len(P) == 3:
        return 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
    return 0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[1]))


def survey(feb_path):
    parts, inst, asm_nodes, conns, behav = parse_inp()
    f = Feb(feb_path)
    ni = NodeIndex(f)
    elem_of = defaultdict(set)
    for name, (et, d) in f.elem_blocks.items():
        for conn in d.values():
            for nd in conn:
                elem_of[nd].add(name)
    bcs = bc_nodes(f)
    fams = defaultdict(list)
    for c in conns:
        fams[family(c['elset'])].append(c)
    res = {}
    for fam, cs in sorted(fams.items()):
        loft = FAMILY_LOFT.get(fam)
        rows = []
        for c in cs:
            xa, pa = end_xyz(c['a'], parts, inst, asm_nodes)
            xb, pb = end_xyz(c['b'], parts, inst, asm_nodes)
            fa, da = ni.nearest(xa)
            fb, db = ni.nearest(xb)
            rows.append(dict(c, xa=xa, xb=xb, pa=pa, pb=pb, fa=fa, fb=fb, da=da, db=db, L=np.linalg.norm(xb - xa)))
        res[fam] = (loft, rows)
    return f, parts, inst, asm_nodes, conns, behav, res, elem_of, bcs


def main():
    feb_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RUNS_DIR, 'LPFmsTk_DM1', 'LPFmsTk_DM1.feb')
    f, parts, inst, asm_nodes, conns, behav, res, elem_of, bcs = survey(feb_path)
    print(f'{len(conns)} CONN3D2 connectors; {sum(1 for c in conns if c.get("elset_mismatch"))} elset mismatches; '
          f'model {os.path.basename(feb_path)}')
    for fam, (loft, rows) in res.items():
        beh = Counter(r['behavior'] for r in rows)
        Ls = np.array([r['L'] for r in rows])
        ends_a = Counter(r['pa'] for r in rows)
        print(f'\n== {fam}: {len(rows)} connectors, {dict(ends_a)} -> {dict(Counter(r["pb"] for r in rows))}; '
              f'L {Ls.min():.1f}-{Ls.max():.1f} mm; behaviours {dict(beh)}')
        da = max(r['da'] for r in rows)
        db = max(r['db'] for r in rows)
        uniq_a = len({r['fa'] for r in rows})
        uniq_b = len({r['fb'] for r in rows})
        print(f'   FEBio match: max dist end A {da:.2e}, end B {db:.2e} mm; distinct FEBio nodes A {uniq_a}, B {uniq_b}')
        if not loft:
            print('   (no loft)')
            continue
        d = loft_elements(f, loft)
        lnodes = {n for conn in d.values() for n in conn}
        area = sum(tri_area(f, conn) for conn in d.values())
        ina = sum(r['fa'] in lnodes for r in rows)
        inb = sum(r['fb'] in lnodes for r in rows)
        print(f'   loft {loft}: {len(d)} elements, {len(lnodes)} nodes, area {area:.0f} mm^2; '
              f'connector ends that are loft nodes: A {ina}/{len(rows)}, B {inb}/{len(rows)}')
        # boundary nodes of the loft and what holds them
        ecount = Counter()
        for conn in d.values():
            for k in range(len(conn)):
                ecount[tuple(sorted((conn[k], conn[(k + 1) % len(conn)])))] += 1
        bnodes = {n for e, c in ecount.items() if c == 1 for n in e}
        held = Counter()
        for nd in bnodes:
            others = sorted(elem_of[nd] - {loft})
            tag = ('BC' if nd in bcs else '') + ('+' if nd in bcs and others else '') + (','.join(others[:2]) if others else '')
            held[tag or 'FREE'] += 1
        print(f'   loft boundary nodes {len(bnodes)}; held by: ' + '; '.join(f'{k} {v}' for k, v in held.most_common(8)))
        a_nodes = {r['fa'] for r in rows}
        a_hold = Counter(('BC:' + '/'.join(sorted(set(bcs[n])))[:40] if n in bcs else 'noBC') + ' | ' +
                         ','.join(sorted(elem_of[n] - {loft})[:3]) for n in a_nodes)
        print(f'   anchor (A) ends held by: ' + '; '.join(f'{k} x{v}' for k, v in a_hold.most_common(5)))


if __name__ == '__main__':
    main()
