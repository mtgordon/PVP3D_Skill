"""Batch 234 (2026-10-09, the user: "try running those again with a stiffer cervix"). Two tube runs blew up at the cervix
(all fastest nodes in _PickedSet346): L202_distpress_apexclosed_healthy (t 0.806) and L232_curvepress_apexclosed_pmc_P1 (t 0.873).
On the tube line the walls are Yeoh x0.35 of the source, the cervix with them. The user's rule (2026-10-09): the cervix starts
at 3x the vaginal wall (up to ~10x). This sets Vagina_Cervix's c1, c2, k to FACTOR x the run's own Vagina_AVW values (the
wall), changing nothing else.
usage: py -3.10 build_batch234.py SRC_RUN NEW_NAME FACTOR"""
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
    get = lambda nm: re.search(r'<material id="\d+" name="%s" type="Yeoh">.*?</material>' % nm, t, re.S).group(0)  # noqa: E731
    wall, cx = get('Vagina_AVW'), get('Vagina_Cervix')
    new, log = cx, []
    for k in ('c1', 'c2', 'k'):
        w = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), wall).group(1))
        old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), cx).group(1))
        new = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{w * F:.7g}</{k}>', new)
        log.append(f'{k} {old:.7g} -> {w * F:.7g}')
    t = t.replace(cx, new)
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch234.py)\nNOT FAITHFUL (the user, 2026-10-09: a '
                   f'stiffer cervix): Vagina_Cervix Yeoh = {F:g} x the wall (Vagina_AVW) of this run: ' + ', '.join(log) + '\n')
    print(name, ', '.join(log))
