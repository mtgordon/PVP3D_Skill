"""Batches 172 / 173 (2026-10-05, the user: "stop them and run the tests"): why the refined-curve runs (L171) crawl.
L171_wrapref_endspr_selfall_walls035_P1 ran like L169 to t 0.30, then needed 181 steps for t 0.30-0.35 (L169: 24) with
initial residuals ~300x L169's; from t 0.30 the curves' self-contact engaged more and more faces (2 -> 12 by t 0.36) at rest z
-28 .. -22 mm, just above the refined zone (the 0.67 mm elements meet the ~2 mm ones there). Two single-change tests of P1:
  L172_ref_selfpen1_P1          L171_..._P1 with the self-contact (curve_all_self) penalty 5 -> 1 (auto_penalty kept)
  L173_ref_compose_P1           the refined curves rebuilt with their intermediate rings' end nodes following the column
                                nodes directly (build_batch135.REFINE_COMPOSE: no constraint refers to a constrained node);
                                otherwise as L171 (build_batch171's chain), P1 only
usage: py -3.10 build_batch172.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
import build_batch135 as b135  # noqa: E402
import build_batch166 as b166  # noqa: E402
import build_batch171 as b171  # noqa: E402
from build_batch162 import write  # noqa: E402
from build_batch165 import wall_damping  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L171_wrapref_endspr_selfall_walls035_P1'

if __name__ == '__main__':
    # A: the self-contact penalty
    p = os.path.join(RUNS, SRC, SRC + '.feb')
    t = open(p, encoding='utf-8').read()
    a = t.index('<contact name="curve_all_self"')
    b = t.index('</contact>', a)
    con, k = re.subn(r'<penalty>5</penalty>', '<penalty>1</penalty>', t[a:b])
    assert k == 1
    write('L172_ref_selfpen1_P1', t[:a] + con + t[b:], f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch172.py)\n'
                                                       f'test: the curves\' self-contact penalty 5 -> 1 (auto_penalty kept)\n')
    # B: composed constraints on the refined rings
    b135.REFINE_DISTAL, b135.REFINE_COMPOSE = b171.REFINE, True
    b166.GEOM, b166.OUTER, b166.ENDSPR = 'L173_thin275_wrapend_full_refc', 'L173_wrapref_outer', 'L173_wrapref_outer_endspr'
    b166.step1()
    with open(os.path.join(RUNS, b166.GEOM, b166.GEOM + '.feb.changes.txt'), 'a', encoding='utf-8') as fo:
        fo.write(f'NOT IN SOURCE: the curves\' last {b171.REFINE[0]} station gaps at the distal end split into {b171.REFINE[1]}; '
                 f'test: the intermediate rings\' end nodes follow the column nodes directly (REFINE_COMPOSE; '
                 f'tools/build_batch172.py)\n')
    b166.step2()
    b166.step34()
    q = os.path.join(RUNS, b166.ENDSPR, b166.ENDSPR + '.feb')
    t, msg5 = b171.self_all(open(q, encoding='utf-8').read())
    t, msg6 = wall_damping(t)
    name = 'L173_wrapref_endspr_selfall_damp'
    write(name, t, f'base: claude_diag/runs/{b166.ENDSPR}/{b166.ENDSPR}.feb (tools/build_batch172.py)\n{msg5}\n{msg6}\n')
    print(' '.join(build(173, 'ref_compose', {'walls': 0.35}, ['P1'], name)))
