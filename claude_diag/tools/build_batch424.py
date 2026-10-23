"""Batch 424 (2026-10-22; the user, away 8-10 h: "if this does solve the issue, start trying to get a max and baseline versions").
build_batch244's case switch from ANY case (not only P1): the file's support values are checked against case FROM, then
rewritten to case TO (build_batch244.CASES: levator % of the healthy law, CL / USL scale, anterior AVW-Para scale, posterior
P-arcus scale); nothing else changes. baseline = 'healthy' (all x1, levator 100 %).
usage: py -3.10 build_batch424.py SRC_RUN FROM_CASE NEW_RUN TO_CASE"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch244 import CASES  # noqa: E402
from build_batch107 import LA100  # noqa: E402
from build_batch162 import write  # noqa: E402
from junction_view_data import load  # noqa: E402


def convert(t, frm, to):
    A, B = CASES[frm], CASES[to]
    log = []

    def scale(m):
        nm, old = m.group(2), float(m.group(3))
        if re.match(r'(CL|USL)-[LR]_conn_\d_mat$', nm):
            want, new = A[1], B[1]
        elif re.match(r'AVW-Para-[LR]_conn_mat$', nm):
            want, new = A[2], B[2]
        elif nm == 'Parcus_conn_mat':
            want, new = A[3], B[3]
        else:
            return m.group(0)
        assert abs(old - want) < 1e-9, (nm, old, want)
        log.append(f'{nm} {old:g} -> {new:g}')
        return f'{m.group(1)}{new:g}{m.group(4)}'
    t = re.sub(r'(<discrete_material id="\d+" name="([^"]+)" type="[^"]+">\s*<scale>)([^<]+)(</scale>)', scale, t)
    assert len(log) == 13, log
    m = re.search(r'<material id="\d+" name="LA_Yamada50pct_Yeoh" type="Yeoh">.*?</material>', t, re.S)
    blk = m.group(0)
    for k, v in LA100.items():
        old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), blk).group(1))
        assert abs(old - v * A[0] / 100) < 1e-6 * max(1, v), (k, old)
        blk = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{v * B[0] / 100:.7g}</{k}>', blk)
    t = t.replace(m.group(0), blk)
    log.append(f'levator {A[0]} % -> {B[0]} % of the healthy law')
    return t, log


if __name__ == '__main__':
    src, frm, new, to = sys.argv[1:5]
    t, _, _ = load(src)
    t, log = convert(t, frm, to)
    print(new); print('; '.join(log))
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch424.py)\ncase {frm} -> {to}: ' + '; '.join(log) + '\n')
