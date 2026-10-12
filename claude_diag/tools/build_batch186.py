"""Batch 186 (2026-10-06, the user: "do the fold test first. If the fold isn't right, it will affect the other things").
The tube curves' distal ends fold through themselves from t ~0.8 (curve_all_self, node-on-facet, search_radius 0.5 mm). On the
refined distal mesh (L171), penalty 5 crawled from t 0.30; penalty 1 (L172_ref_selfpen1_P1) ran to t 1 but let the fold through
(35-70 of 944 outer nodes inside at t 0.8-1; unrefined L169 at penalty 5: 2-26). Two single-change tests from L171 / L172, P1:
  L186_ref_selfpen25_P1        L171_..._P1 with the curve_all_self penalty 5 -> 2.5 (between the crawl and the leak)
  L186_ref_selfpen1_2pass_P1   L172_ref_selfpen1_P1 with two_pass 0 -> 1 (each side of the fold also checked as the primary)
Judge with tools/curve_self_pen.py RUN 0.8,0.85,0.9,1 and the pace (steps for t 0.30-0.36: L169 24, L172 27, L171 181).
usage: py -3.10 build_batch186.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')


def edit(src, tag, new):
    p = os.path.join(RUNS, src, src + '.feb')
    t = open(p, encoding='utf-8').read()
    a = t.index('<contact name="curve_all_self"')
    b = t.index('</contact>', a)
    con, k = re.subn(rf'<{tag}>[^<]*</{tag}>', f'<{tag}>{new}</{tag}>', t[a:b])
    assert k == 1
    old = re.search(rf'<{tag}>([^<]*)</{tag}>', t[a:b]).group(1)
    return t[:a] + con + t[b:], f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch186.py)\ntest: curve_all_self {tag} {old} -> {new}\n'


if __name__ == '__main__':
    write('L186_ref_selfpen25_P1', *edit('L171_wrapref_endspr_selfall_walls035_P1', 'penalty', '2.5'))
    write('L186_ref_selfpen1_2pass_P1', *edit('L172_ref_selfpen1_P1', 'two_pass', '1'))
