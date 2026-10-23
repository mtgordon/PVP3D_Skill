"""Batch 326 (2026-10-18; the user: "is it reasonable that anytime there is jitters to add some mass/damping?"; memory
feedback-jitter-mass-damping). NOT IN SOURCE: the domains DOMS get density x DENS on their materials (copied if another domain
shares the material) and mass damping C CD (body_load, elem_set = the domain). For a crawl traced to a few chattering nodes
(e.g. L312b from t 0.637: AVW-Para solid-loft nodes on the PM's free edge).
usage: py -3.10 build_batch326.py SRC_RUN NEW_RUN --doms A,B [--dens 100] [--cd 200]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new, doms, dens, cd):
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    dm = {d: re.search(rf'<(?:Shell|Solid)Domain name="{re.escape(d)}" mat="([^"]+)"', t).group(1) for d in doms}
    users = {m: re.findall(rf'<(?:Shell|Solid)Domain name="([^"]+)" mat="{re.escape(m)}"', t) for m in set(dm.values())}
    mid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
    done = {}
    for d, m in dm.items():
        if m in done:
            new_m = done[m]
        else:
            mm = re.search(rf'\t\t<material id="\d+" name="{re.escape(m)}" type="[^"]+">.*?</material>\n', t, re.S).group(0)
            heavy = re.sub(r'<density>([^<]+)</density>', lambda q: f'<density>{float(q.group(1)) * dens:.6g}</density>', mm)
            assert heavy != mm, m
            if set(users[m]) <= set(doms):
                t = t.replace(mm, heavy, 1); new_m = m
            else:
                new_m = m + '_heavy'
                heavy = re.sub(r'<material id="\d+" name="[^"]+"', f'<material id="{mid}" name="{new_m}"', heavy); mid += 1
                i = t.index('</Material>'); t = t[:i] + heavy + t[i:]
            done[m] = new_m
        if new_m != m:
            t = re.sub(rf'(<(?:Shell|Solid)Domain name="{re.escape(d)}" mat=")[^"]+"', rf'\g<1>{new_m}"', t, count=1)
        i = t.index('</Loads>')
        t = t[:i] + f'\t\t<body_load name="jitter_damping_{d}" type="mass damping" elem_set="{d}">\n\t\t\t<C>{cd:g}</C>\n\t\t</body_load>\n' + t[i:]
    msg = f'NOT IN SOURCE (convergence aid): {", ".join(doms)}: density x{dens:g} (materials {sorted(set(done.values()))}) + mass damping C {cd:g}'
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch326.py)\n{msg}\n')


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda k, d: float(a[a.index(k) + 1]) if k in a else d   # noqa: E731
    main(a[0], a[1], a[a.index('--doms') + 1].split(','), g('--dens', 100.0), g('--cd', 200.0))
