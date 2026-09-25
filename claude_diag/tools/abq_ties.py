"""Resolve the Abaqus *Tie constraints between LA-new and the ATLA / Posterior_Arcus trusses.

Abaqus: slave = LA-new nodes (node set or element-surface nodes), master = truss nodes
(node-based surface), position tolerance 4 mm, adjust=no  ->  each slave node within
tolerance is tied (translations) to its closest master node.
Returns FEBio node ids (by exact coordinate match against LA element nodes / truss nodes).
"""
import re
import numpy as np

from inpmap import part_nodes, INP

TOL = 4.0


def _lines():
    return open(INP, encoding='latin-1').read().splitlines()


def nset(name, kind='nset'):
    out = []
    grab = gen = False
    for ln in _lines():
        s = ln.strip()
        if s.lower().startswith('*' + kind) and re.search(r'%s=%s\s*(,|$)' % (kind, re.escape(name)), s):
            grab = True
            gen = 'generate' in s.lower()
            continue
        if grab:
            if s.startswith('*'):
                break
            v = [int(x) for x in s.replace(' ', '').split(',') if x]
            out += list(range(v[0], v[1] + 1, v[2] if len(v) > 2 else 1)) if gen else v
    return out


def la_elements():
    """LA-new element connectivity {eid: [local node ids]}."""
    out = {}
    inpart = inel = False
    for ln in _lines():
        s = ln.strip()
        if s.lower().startswith('*part,'):
            inpart = 'name=LA-new' in s
            continue
        if not inpart:
            continue
        if s.lower().startswith('*end part'):
            break
        if s.startswith('*'):
            inel = s.lower().startswith('*element')
            continue
        if inel and s:
            v = [int(x) for x in s.replace(' ', '').split(',') if x]
            out[v[0]] = v[1:]
    return out


TIES = {
    # name: (slave spec, master part)
    'ATLA_L': (('nset', '_PickedSet473'), 'ATLA_Left'),
    'ATLA_R': (('nset', '_PickedSet475'), 'ATLA_Right'),
    'PArcus_L': (('elset', '__PickedSurf737_SPOS'), 'Posterior_Arcus _Left'),
    'PArcus_R': (('elset', '__PickedSurf739_SPOS'), 'Posterior_Arcus _Right'),
}


def resolve():
    la = part_nodes('LA-new')
    els = la_elements()
    res = {}
    for name, ((kind, setname), mpart) in TIES.items():
        if kind == 'nset':
            slaves = nset(setname, 'nset')
        else:
            slaves = sorted({n for e in nset(setname, 'elset') for n in els[e]})
        mn = part_nodes(mpart)
        mids = sorted(mn)
        M = np.array([mn[i] for i in mids])
        pairs = []
        for s in slaves:
            d = np.linalg.norm(M - la[s], axis=1)
            k = int(np.argmin(d))
            if d[k] <= TOL:
                pairs.append((s, mids[k], float(d[k])))
        res[name] = (pairs, mn)
    return res, la


if __name__ == '__main__':
    res, la = resolve()
    for name, (pairs, mn) in res.items():
        used = sorted({m for _, m, _ in pairs})
        free = [m for m in sorted(mn) if m not in used]
        ds = [d for *_, d in pairs]
        print(f'{name}: {len(pairs)} LA slave nodes tied to {len(used)} of {len(mn)} truss nodes; '
              f'dist min/med/max {min(ds):.3f}/{np.median(ds):.3f}/{max(ds):.3f}; truss nodes with no slave: {free}')
