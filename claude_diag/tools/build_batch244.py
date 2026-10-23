"""Batch 244 (2026-10-09, the user's range-of-motion goal on BOTH lines; the tube line's models exist only as P1 files). As
build_batch231, but for every build_batch243 case: the P1 file's case values (levator 80 % of the healthy law, every CL / USL
<scale> 0.7, Parcus 0.15, anterior AVW-Para 1) are checked, then rewritten to the case's values; nothing else changes.
  case     levator %  apical (CL / USL)  anterior (AVW-Para)  posterior (Parcus)
  healthy  100        1                  1                    1
  P1       80         0.7                1                    0.15
  P2       40         0.4                1                    0.15
  F7A      50         0.2                0.5                  0.1      (Fig. 7A rectocele)
  F8A      50         1                  0.1                  0.5      (Fig. 8A cystocele)
  F6a0     80         1                  1                    0.15     (Fig. 6, no apical impairment)
  F6a60    80         0.4                1                    0.15     (Fig. 6, apical 60 %)
usage: py -3.10 build_batch244.py SRC_P1_RUN NEW_PREFIX CASE[,CASE...]   -> runs/<NEW_PREFIX>_<CASE>/"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch107 import LA100  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
CASES = {'healthy': (100, 1, 1, 1), 'P1': (80, 0.7, 1, 0.15), 'P2': (40, 0.4, 1, 0.15), 'F7A': (50, 0.2, 0.5, 0.1),
         'F8A': (50, 1, 0.1, 0.5), 'F6a0': (80, 1, 1, 0.15), 'F6a60': (80, 0.4, 1, 0.15),
         'MAX': (40, 0.2, 0.5, 0.05),   # 2026-10-09 evening: the heaviest combination (levator 60 %, apical 80 %, anterior 50 %, posterior 95 %)
         # 2026-10-18 MAX-lite (the user: MAX goes beyond what is needed): MAX with ONE support less impaired each
         'MAXl_post02': (40, 0.2, 0.5, 0.2), 'MAXl_ap05': (40, 0.5, 0.5, 0.05), 'MAXl_la60': (60, 0.2, 0.5, 0.05),
         'MAXl_ant1': (40, 0.2, 1, 0.05),
         # per-target / second level: posterior unimpaired (anterior-heavy, for Ba), levator 50 %
         'MAXl_post1': (40, 0.2, 0.5, 1), 'MAXl_la50': (50, 0.2, 0.5, 0.05),
         'MAXl_post02ap05': (40, 0.5, 0.5, 0.2),   # the two that delay the snap, together
         'MAXl_post01': (40, 0.2, 0.5, 0.1), 'MAXl_post0075': (40, 0.2, 0.5, 0.075)}   # inside the P-arcus cliff (x0.05 prolapses, x0.2 does not)
P1 = CASES['P1']


def convert(t, case):
    la, ap, ant, post = CASES[case]
    log = []

    def scale(m):
        nm, old = m.group(2), float(m.group(3))
        if re.match(r'(CL|USL)-[LR]_conn_\d_mat$', nm):
            want, new = P1[1], ap
        elif re.match(r'AVW-Para-[LR]_conn_mat$', nm):
            want, new = P1[2], ant
        elif nm == 'Parcus_conn_mat':
            want, new = P1[3], post
        else:
            return m.group(0)
        assert abs(old - want) < 1e-9, (nm, old)
        log.append(f'{nm} scale {old:g} -> {new:g}')
        return f'{m.group(1)}{new:g}{m.group(4)}'
    t = re.sub(r'(<discrete_material id="\d+" name="([^"]+)" type="[^"]+">\s*<scale>)([^<]+)(</scale>)', scale, t)
    assert len(log) == 13, log
    m = re.search(r'<material id="\d+" name="LA_Yamada50pct_Yeoh" type="Yeoh">.*?</material>', t, re.S)
    blk = m.group(0)
    for k, v in LA100.items():
        old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), blk).group(1))
        assert abs(old - v * P1[0] / 100) < 1e-6 * max(1, v), (k, old)
        blk = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{v * la / 100:.7g}</{k}>', blk)
    t = t.replace(m.group(0), blk)
    log.append(f'levator LA_Yamada50pct_Yeoh c1, c2, k at {la} % of the healthy law (was {P1[0]} %)')
    return t, log


if __name__ == '__main__':
    src, prefix, cases = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
    t0 = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    for c in cases:
        t, log = convert(t0, c)
        write(f'{prefix}_{c}', t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch244.py)\nthe paper\'s case {c} '
                                  f'(Luo et al. 2015, the range-of-motion checks) instead of P1: ' + '; '.join(log) + '\n')
