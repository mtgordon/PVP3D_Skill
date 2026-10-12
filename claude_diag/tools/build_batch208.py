"""Batch 208 (2026-10-07 evening). The user: the seam lips "should affect the run as little as possible as they are not in the
original design". Check how much they do: every lip (AVW_ext, cervix_ext, PVW_ext) on its own copy of its wall's Yeoh with
c1, c2, k x F (as build_batch145.wall_scale), every original wall element unchanged. If the answer moves little from the
lips at the wall's own material, the lips are neutral.
usage: py -3.10 build_batch208.py SRC_RUN NEW_NAME FACTOR"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
LIPS = {'AVW_ext': 'Vagina_AVW', 'cervix_ext': 'Vagina_Cervix', 'PVW_ext': 'Vagina_PVW'}

if __name__ == '__main__':
    src, name, F = sys.argv[1], sys.argv[2], float(sys.argv[3])
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    a = t.index('<Material>')
    b = t.index('</Material>', a)
    nid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t[a:b])) + 1
    add, log = '', []
    for dn, wall in LIPS.items():
        m = re.search(r'<material id="\d+" name="%s" type="Yeoh">.*?</material>' % re.escape(wall), t, re.S).group(0)
        new = re.sub(r'<material id="\d+" name="[^"]+"', f'<material id="{nid}" name="{dn}_lip"', m)
        for k in ('c1', 'c2', 'k'):
            v = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), m).group(1))
            new = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{v * F:.7g}</{k}>', new)
        add += '\t' + new + '\n\t'
        t, n = re.subn(r'(<SolidDomain name="%s" mat=")[^"]+(")' % dn, r'\g<1>%s_lip\2' % dn, t)
        assert n == 1, dn
        log.append(f'{dn} -> {dn}_lip')
        nid += 1
    b = t.index('</Material>')
    t = t[:b] + add + t[b:]
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch208.py)\n'
                   f'lip sensitivity (lips NOT IN SOURCE): every lip on its own copy of its wall\'s Yeoh, c1, c2, k x{F:g} '
                   f'({"; ".join(log)}); every original wall element unchanged\n')
