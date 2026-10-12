"""Batch 217 (2026-10-08, unattended; DIAGNOSTIC). The AVW's paravaginal springs (AVW-Para-L/R_conn, 35 a side to one pinned
point) shorten as the AVW descends and, tension-only as in the source, carry nothing (README_2026-10-08). With compression at
the table's first slope (0.58333 N / 4.7 mm = 0.124 N/mm, L215_src_sl01_paracomp_P1) the AVW stays 36 mm higher. How little
compression is enough? Here the tension side is the source table unchanged; one point (-47.115, -F * f) is added before it so
the compression stiffness is f x the first slope (linear, constant beyond -47 mm). NOT IN SOURCE.
usage: py -3.10 build_batch217.py SRC_RUN NEW_NAME FRACTION"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')

if __name__ == '__main__':
    src, name, f = sys.argv[1], sys.argv[2], float(sys.argv[3])
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    for mat in ('AVW-Para-L_conn_mat', 'AVW-Para-R_conn_mat'):
        m = re.search(r'<discrete_material id="\d+" name="%s" type="nonlinear spring">.*?</discrete_material>' % re.escape(mat),
                      t, re.S)
        blk = m.group(0)
        pts = [tuple(float(v) for v in p.split(',')) for p in re.findall(r'<pt>([^<]+)</pt>', blk)]
        assert pts[0] == (0.0, 0.0), pts[0]
        k1 = pts[1][1] / pts[1][0]
        x = -pts[-1][0]
        new = blk.replace('<points>', f'<points>\n\t\t\t\t\t<pt>{x:g},{k1 * f * x:.6g}</pt>', 1)
        assert '<extend>constant</extend>' in new
        t = t.replace(blk, new)
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch217.py)\n'
                   f'DIAGNOSTIC, NOT IN SOURCE: AVW-Para-L/R_conn carry compression at {f:g} x the first slope of their source '
                   f'table ({k1 * f:.4g} N/mm per spring; point ({x:g}, {k1 * f * x:.4g}) added; tension side unchanged)\n')
