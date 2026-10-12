"""Batch 231 (2026-10-09, the user: "run healthy and P2 for L202_distpress_apexclosed"). The tube variant was built only as a
P1 file. Its P1-specific values (build_batch107 case P1 on healthy supports at the source values): levator at 80 % of the
healthy Yamada law (LA_Yamada50pct_Yeoh c1, c2, k = 0.8 x LA100), every CL / USL set's <scale> 0.7, Parcus_conn_mat 0.15.
This rewrites exactly those to the other case and changes nothing else:
  healthy: levator 100 %, CL / USL 1, Parcus 1        P2: levator 40 %, CL / USL 0.4, Parcus 0.15
usage: py -3.10 build_batch231.py SRC_P1_RUN NEW_NAME CASE"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch107 import LA100  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
CASES = {'healthy': (100, 1.0, 1.0), 'P2': (40, 0.4, 0.15)}
P1 = (80, 0.7, 0.15)

if __name__ == '__main__':
    src, name, case = sys.argv[1:4]
    la, clusl, parcus = CASES[case]
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    log = []
    # connectors: CL / USL sets and Parcus, checked to hold the P1 values first
    def scale(m):
        nm, old = m.group(2), float(m.group(3))
        if re.match(r'(CL|USL)-[LR]_conn_\d_mat$', nm):
            assert abs(old - P1[1]) < 1e-9, (nm, old); new = clusl
        elif nm == 'Parcus_conn_mat':
            assert abs(old - P1[2]) < 1e-9, (nm, old); new = parcus
        else:
            return m.group(0)
        log.append(f'{nm} scale {old:g} -> {new:g}')
        return f'{m.group(1)}{new:g}{m.group(4)}'
    t = re.sub(r'(<discrete_material id="\d+" name="([^"]+)" type="[^"]+">\s*<scale>)([^<]+)(</scale>)', scale, t)
    assert len(log) == 11, log
    # levator
    m = re.search(r'<material id="\d+" name="LA_Yamada50pct_Yeoh" type="Yeoh">.*?</material>', t, re.S)
    blk = m.group(0)
    for k, v in LA100.items():
        old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), blk).group(1))
        assert abs(old - v * P1[0] / 100) < 1e-6 * max(1, v), (k, old)
        blk = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{v * la / 100:.7g}</{k}>', blk)
    t = t.replace(m.group(0), blk)
    log.append(f'levator LA_Yamada50pct_Yeoh c1, c2, k at {la} % of the healthy law (was {P1[0]} %)')
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch231.py)\nthe paper\'s case {case} instead of P1 '
                   f'(Luo et al. 2015; build_batch107 values): ' + '; '.join(log) + '\n')
    print(case, '; '.join(log[-2:]))
