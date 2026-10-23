"""Batch 323 (2026-10-18 night). L309 (every loft x0.01 + density x100 + loft mass damping, with L298's edge spring patches)
does NOT snap at t 0.73 the way every chains run with full lofts does (L294b / L298 / L310 / L314 / L317); its folds creep
instead (0.69-0.93 mm at t 0.76), like the springs line with patches (L288c). So the loft sheets' pull on single curve nodes
looks like what triggers the snap. L306 (those triangles removed) crawled. Here the local form of L309 (NOT IN SOURCE):
every loft-sheet tri3 with a node on the distal curves (curve nodes at rest z < ZCUT) moves to its own shell domain
<sheet>_cb with its own material (the sheet's, E x SCALE, density x DENS) and mass damping C CD; the rest of every loft is
unchanged (full stiffness, no extra mass).
usage: py -3.10 build_batch323.py SRC_RUN NEW_RUN [--zcut -15] [--scale 0.01] [--dens 100] [--cd 200]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, zcut, scale, dens, cd):
    t, X, E = load(src)
    curve = {v for nm, c in E.values() if nm.startswith('canal_wrap_hex') and nm != 'canal_wrap_hex_avw_cervix' for v in c}
    dist = {n for n in curve if X[n][2] < zcut}
    sheets = sorted({nm for nm, c in E.values() if nm.endswith('_sheet') and len(c) == 3})
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    counts = {}
    for sh in sheets:
        ids = [k for k, (nm, c) in E.items() if nm == sh and any(v in dist for v in c)]
        if not ids:
            continue
        a = t.index(f'<Elements type="tri3" name="{sh}">'); b = t.index('</Elements>', a)
        blk = t[a:b]; rows = []
        for k in ids:
            row = f'\t\t\t<elem id="{k}">' + ','.join(map(str, E[k][1])) + '</elem>\n'
            assert row in blk, (sh, k)
            blk = blk.replace(row, '', 1); rows.append(row)
        assert '<elem ' in blk, f'{sh}: every triangle touches the curve end'
        b2 = b + len('</Elements>\n')
        t = t[:a] + blk + t[b:b2] + f'\t\t<Elements type="tri3" name="{sh}_cb">\n' + ''.join(rows) + '\t\t</Elements>\n' + t[b2:]
        # domain + material
        md = re.search(rf'\t\t<ShellDomain name="{sh}" mat="([^"]+)">.*?</ShellDomain>\n', t, re.S)
        mat = md.group(1)
        t = t.replace(md.group(0), md.group(0) + md.group(0).replace(f'name="{sh}" mat="{mat}"', f'name="{sh}_cb" mat="{mat}_cb"'), 1)
        mm = re.search(rf'\t\t<material id="\d+" name="{re.escape(mat)}" type="neo-Hookean">.*?</material>\n', t, re.S)
        m2 = re.sub(r'<material id="\d+" name="[^"]+"', f'<material id="{mid}" name="{mat}_cb"', mm.group(0)); mid += 1
        m2 = re.sub(r'<E>([^<]+)</E>', lambda m: f'<E>{float(m.group(1)) * scale:.6g}</E>', m2)
        m2 = re.sub(r'<density>([^<]+)</density>', lambda m: f'<density>{float(m.group(1)) * dens:.6g}</density>', m2)
        i = t.index('</Material>')
        t = t[:i] + m2 + t[i:]
        i = t.index('</Loads>')
        t = t[:i] + f'\t\t<body_load name="loft_damping_{sh}_cb" type="mass damping" elem_set="{sh}_cb">\n\t\t\t<C>{cd:g}</C>\n\t\t</body_load>\n' + t[i:]
        counts[sh] = len(rows)
    msg = (f'NOT IN SOURCE: the loft tri3 touching a distal curve node (rest z < {zcut}) moved to their own domains <sheet>_cb {counts} '
           f'with their own materials (E x{scale:g}, density x{dens:g}) and mass damping C {cd:g}; the rest of every loft unchanged')
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch323.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], g('--zcut', -15.0), g('--scale', 0.01), g('--dens', 100.0), g('--cd', 200.0))
