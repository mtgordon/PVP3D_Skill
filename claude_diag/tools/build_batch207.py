"""Batch 207 (2026-10-07 evening, the user: "is there anything you can do to get the lips line to give the paper data?").
On the lips line (L205_ext5_*: 5 mm unloaded lips past the contact-only seam, no band) the AVW drops without the band (P1 gap
Bp - Ba -3.3 vs the band default's +3.9). On the band line the gap came from a stiff strip along the AVW / cervix side of the
seam with the PVW strip-free. Here the strip goes into the ADDED lip elements only: the AVW + cervix lips (AVW_ext,
cervix_ext) get their own Yeoh, c1 x F (c2, k, density as the wall, as the band scaled c1 only); the PVW lips and every
original wall element keep the source's material. NOT IN SOURCE (local, added material only).
usage: py -3.10 build_batch207.py SRC_RUN NEW_NAME FACTOR"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')

if __name__ == '__main__':
    src, name, F = sys.argv[1], sys.argv[2], float(sys.argv[3])
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    m = re.search(r'<material id="\d+" name="Vagina_AVW" type="Yeoh">.*?</material>', t, re.S).group(0)
    c1 = float(re.search(r'<c1>([^<]+)</c1>', m).group(1))
    a = t.index('<Material>')
    b = t.index('</Material>', a)
    ids = [int(v) for v in re.findall(r'<material id="(\d+)"', t[a:b])]
    nid = max(ids) + 1
    assert nid == len(ids) + 1, 'material ids not consecutive'
    lip = re.sub(r'<material id="\d+" name="Vagina_AVW"', f'<material id="{nid}" name="AVW_cervix_lip"', m)
    lip = lip.replace(f'<c1>{re.search(r"<c1>([^<]+)</c1>", m).group(1)}</c1>', f'<c1>{c1 * F:.7g}</c1>')
    t = t[:b] + '\t' + lip + '\n\t' + t[b:]
    n = 0
    for dn in ('AVW_ext', 'cervix_ext'):
        t, k = re.subn(r'(<SolidDomain name="%s" mat=")[^"]+(")' % dn, r'\1AVW_cervix_lip\2', t)
        n += k
    assert n == 2
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch207.py)\n'
                   f'NOT IN SOURCE (stiff AVW-side lips): the AVW + cervix lips (AVW_ext, cervix_ext) on their own Yeoh '
                   f'AVW_cervix_lip, c1 x{F:g} ({c1 * F:.7g}; c2, k, density as Vagina_AVW); the PVW lips and every original '
                   f'wall element keep the source material\n')
