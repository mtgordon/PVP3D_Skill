"""Batch 338 (2026-10-19; the user, on L326d (C 50) jittering at the loft <-> PM interaction and L326b (C 200) lagging: damping
only on a band along the PM free edge). From a build_batch326 run (domains DOMS with heavy materials + jitter_damping body loads):
the elements of DOMS with a node within R mm (at rest) of the PM's FREE edge (PM_Plane boundary nodes off BC-PM-outer-arc) move
to their own domain <dom>_band (same heavy material, same shell settings) and keep the mass damping (C as in SRC); the rest of
each domain gets a LIGHT copy of its material (density / DENS, i.e. back to normal) and no damping. NOT IN SOURCE (convergence aid).
usage: py -3.10 build_batch338.py SRC_RUN NEW_RUN [--r 4] [--dens 100]"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, r=4.0, dens=100.0):
    t, X, E = load(src)
    doms = re.findall(r'<body_load name="jitter_damping_([^"]+)" type="mass damping" elem_set="\1">', t)
    assert doms, 'no jitter_damping loads in ' + src
    # PM free-edge nodes
    pm = [c for nm, c in E.values() if nm == 'PM_Plane']
    fixed = {int(v) for v in re.findall(r'\d+', re.search(r'<NodeSet name="BC-PM-outer-arc">(.*?)</NodeSet>', t, re.S).group(1))}
    cnt = Counter(frozenset((q[i], q[(i + 1) % len(q)])) for q in pm for i in range(len(q)))
    free = {v for e, k in cnt.items() if k == 1 for v in e} - fixed
    F = np.array([X[v] for v in sorted(free)])
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    log = [f'PM free-edge nodes: {len(free)} (PM_Plane boundary off BC-PM-outer-arc)']
    for d in doms:
        m = re.search(rf'(\t\t<Elements type="([^"]+)" name="{re.escape(d)}">\n)(.*?)(\t\t</Elements>\n)', t, re.S)
        rows = re.findall(r'\t\t\t<elem id="\d+">[0-9,]+</elem>\n', m.group(3))
        band, rest = [], []
        for row in rows:
            ns = [int(v) for v in re.search(r'>([0-9,]+)<', row).group(1).split(',')]
            dmin = min(np.min(np.linalg.norm(F - X[v], axis=1)) for v in ns)
            (band if dmin <= r else rest).append(row)
        log.append(f'{d}: {len(band)} of {len(rows)} elements in the band')
        if band and not rest:   # (L338c fix) the whole domain is in the band: leave it as it is (heavy, damped)
            log[-1] += ' (all: left heavy and damped, unchanged)'
            continue
        if not band:   # nothing near the edge: the whole domain goes light, no damping
            rest_rows = rows
        else:
            rest_rows = rest
        dom = re.search(rf'\t\t<(Shell|Solid)Domain name="{re.escape(d)}" mat="([^"]+)"(.*?)(/>\n|>\n.*?</\1Domain>\n)', t, re.S)
        heavy = dom.group(2)
        mm = re.search(rf'\t\t<material id="\d+" name="{re.escape(heavy)}" type="[^"]+">.*?</material>\n', t, re.S).group(0)
        light_name = heavy + '_light'
        if light_name not in t:
            light = re.sub(r'<density>([^<]+)</density>', lambda q: f'<density>{float(q.group(1)) / dens:.6g}</density>', mm)
            light = re.sub(r'<material id="\d+" name="[^"]+"', f'<material id="{mid}" name="{light_name}"', light); mid += 1
            i = t.index('</Material>'); t = t[:i] + light + t[i:]
        # elements: rest stays in d (light), band to d_band (heavy)
        newblk = m.group(1) + ''.join(rest_rows) + m.group(4)
        if band and rest_rows:
            newblk += m.group(1).replace(f'name="{d}"', f'name="{d}_band"') + ''.join(band) + m.group(4)
        t = t.replace(m.group(0), newblk, 1)
        # (2026-10-19 night, L343: fibre lofts) per-element data on d (e.g. the fibre mat_axis) is listed by position in the set:
        # split it as the elements were split
        for ed in re.findall(rf'[ \t]*<ElementData type="[^"]+" elem_set="{re.escape(d)}">\n.*?</ElementData>\n', t, re.S):
            erows = re.findall(r'\t\t\t<elem lid="\d+">(.*?)</elem>\n', ed, re.S)
            assert len(erows) == len(rows), (d, len(erows), len(rows))
            inband = {id(rw) for rw in band}
            er_rest = [e_ for rw, e_ in zip(rows, erows) if id(rw) not in inband]
            er_band = [e_ for rw, e_ in zip(rows, erows) if id(rw) in inband]
            head = re.match(r'[ \t]*<ElementData[^>]*>\n', ed).group(0).lstrip()
            mk = lambda name, er: ('\t\t' + head.replace(f'elem_set="{d}"', f'elem_set="{name}"') +   # noqa: E731
                                   ''.join(f'\t\t\t<elem lid="{k}">{e_}</elem>\n' for k, e_ in enumerate(er, 1)) + '\t\t</ElementData>\n')
            t = t.replace(ed, mk(d, er_rest) + (mk(f'{d}_band', er_band) if band and rest_rows else ''), 1)
            log.append(f'{d}: its {head.split(chr(34))[1]} ElementData split {len(er_rest)} / {len(er_band)}')
        dom_old = dom.group(0)
        dom_light = dom_old.replace(f'mat="{heavy}"', f'mat="{light_name}"', 1)
        dom_band = dom_old.replace(f'name="{d}"', f'name="{d}_band"', 1) if band and rest_rows else ''
        t = t.replace(dom_old, dom_light + dom_band, 1)
        bl = re.search(rf'\t\t<body_load name="jitter_damping_{re.escape(d)}" type="mass damping" elem_set="{re.escape(d)}">.*?</body_load>\n', t, re.S).group(0)
        if band and rest_rows:
            t = t.replace(bl, bl.replace(f'name="jitter_damping_{d}"', f'name="jitter_damping_{d}_band"').replace(f'elem_set="{d}"', f'elem_set="{d}_band"'), 1)
        elif band:   # the whole domain is in the band: keep it heavy and damped as it was
            t = t.replace(dom_light, dom_old, 1)
        else:
            t = t.replace(bl, '', 1)
    msg = (f'NOT IN SOURCE (convergence aid): the jitter mass / damping of {src} kept only on a band within {r:g} mm (at rest) of the '
           f'PM free edge (own domains *_band, heavy material, mass damping as before); the rest of each domain back to normal '
           f'density (*_light materials, density / {dens:g}) and undamped. ' + '; '.join(log))
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch338.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--r', 4.0), g('--dens', 100.0))
