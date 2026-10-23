"""Batch 306 (2026-10-17; the user: "I'm not sure why the chains with lofts fail when the just chains don't. Can't we turn the
lofts stiffness down so much that it basically has no effect?"). Turning the lofts down x0.1 / x0.01 (L305_*) made the chain
nodes nearly free (they are held sideways only by the loft sheet) and the runs crawled from the start. So instead (NOT IN
SOURCE): every loft-sheet tri3 with a node on the distal curves (curve nodes at rest z < ZCUT) is removed, so the curve end
feels only the chain springs (patched as in L288c, when the source has them) and not the sheets' pull on single curve nodes;
the rest of each sheet still holds its chain nodes. Reports any chain node left in no sheet triangle (would be near-free).
usage: py -3.10 build_batch306.py SRC_RUN NEW_RUN [--zcut -15]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, zcut):
    t, X, E = load(src)
    curve = {v for nm, c in E.values() if nm.startswith('canal_wrap_hex') and nm != 'canal_wrap_hex_avw_cervix' for v in c}
    dist = {n for n in curve if X[n][2] < zcut}
    sheets = sorted({nm for nm, c in E.values() if nm.endswith('_sheet') and len(c) == 3})
    cut, note = set(), []
    for s in sheets:
        k = [e for e, (nm, c) in E.items() if nm == s and any(v in dist for v in c)]
        cut |= set(k)
        if k:
            note.append(f'{s} {len(k)} of {sum(1 for nm, _ in E.values() if nm == s)}')
    keep_nodes = {v for e, (nm, c) in E.items() if nm in sheets and e not in cut for v in c}
    chain = set()
    for blk in re.findall(r'<Nodes name="[^"]*_chain_nodes">(.*?)</Nodes>', t, re.S):
        chain |= {int(v) for v in re.findall(r'<node id="(\d+)"', blk)}
    orphan = sorted(n for n in chain if n not in keep_nodes)
    for e in cut:
        row = re.search(rf'\t\t\t<elem id="{e}">[^<]*</elem>\n', t)
        t = t.replace(row.group(0), '', 1)
    msg = (f'NOT IN SOURCE: {len(cut)} loft-sheet tri3 with a node on the distal curves (rest z < {zcut}) removed ({"; ".join(note)}); '
           f'the curve end feels only the chain springs; chain nodes left in no sheet triangle: {len(orphan)}')
    print(new); print(msg)
    if orphan:
        print('orphan chain nodes:', orphan[:40])
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch306.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[a.index('--zcut') + 1]) if '--zcut' in a else -15.0)
