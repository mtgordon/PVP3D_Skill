"""Batch 171 (2026-10-05, the user: "refine the curves at the distal end").
In L169_wrapfull_endspr_selfall_walls035_P2 the curves' distal ends fold back over their own neighbouring elements (within one
or two elements, 1.6-2.8 mm long) from t ~0.75 and then work their way into the curve element by element (nodes cross
internal faces; the steps are already ~3e-4, so smaller steps cannot help): node-on-facet self-contact cannot resolve a
fold that tight at this mesh size. So the curves are rebuilt with the last 6 station gaps at the distal (perineal-body) end
split into 3 along the canal (~1.6 -> ~0.55 mm, about the through-thickness layer); around the arc stays 6 segments, so
the distal end ring and its 5 springs a side are as before (build_batch135.REFINE_DISTAL = (6, 3)).
Otherwise as L169: build_batch166 steps 1-4 (the right side's distal seam node, PM_conn on the outer arc, the 4
PM_avw_bottom springs and 6 new curve_end_PM_conn springs between the curves' distal end ring and the PM's outer arc),
then the self-contact on the curves' whole boundary (build_batch169.all_faces, search_radius 0.5), the wall damping, walls
x0.35, P1 / P2.
  L171_thin275_wrapend_full_ref          the geometry (PM_conn inner arc, as L143)
  L171_wrapref_outer, L171_wrapref_outer_endspr
  L171_wrapref_endspr_selfall_damp       (source wall stiffness)
  L171_wrapref_endspr_selfall_walls035_{P1,P2}
usage: py -3.10 build_batch171.py"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
import build_batch135 as b135  # noqa: E402
import build_batch166 as b166  # noqa: E402
from build_batch164 import add_self_contact  # noqa: E402
from build_batch169 import all_faces, SEARCH_RADIUS  # noqa: E402
from build_batch162 import write, block  # noqa: E402
from build_batch165 import wall_damping  # noqa: E402
from build_batch145 import build  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
REFINE = (6, 3)
b166.GEOM = 'L171_thin275_wrapend_full_ref'
b166.OUTER = 'L171_wrapref_outer'
b166.ENDSPR = 'L171_wrapref_outer_endspr'
NAME = 'L171_wrapref_endspr_selfall_damp'


def self_all(t):
    t, msg5 = add_self_contact(t)                    # the contact block and pair (outer faces), then the whole boundary
    faces, cats = all_faces(t)
    i, j = block(t, 'Surface', 'curve_outer')
    t = (t[:i] + '<Surface name="curve_all">\n'
         + ''.join(f'\t\t\t<quad4 id="{q + 1}">{",".join(map(str, f))}</quad4>\n' for q, f in enumerate(faces)) + '\t\t</Surface>'
         + t[j:])
    for a, b in (('<primary>curve_outer</primary>', '<primary>curve_all</primary>'),
                 ('<secondary>curve_outer</secondary>', '<secondary>curve_all</secondary>'),
                 ('<SurfacePair name="curve_outer_self">', '<SurfacePair name="curve_all_self">'),
                 ('name="curve_outer_self" surface_pair="curve_outer_self"', 'name="curve_all_self" surface_pair="curve_all_self"')):
        assert t.count(a) == 1, a
        t = t.replace(a, b)
    a = t.index('<contact name="curve_all_self"')
    b = t.index('</contact>', a)
    con, k = re.subn(r'<search_radius>[^<]*</search_radius>', f'<search_radius>{SEARCH_RADIUS}</search_radius>', t[a:b])
    assert k == 1
    return t[:a] + con + t[b:], (f'NOT IN SOURCE: the curves\' self-contact on their whole boundary (surface curve_all, '
                                 f'{len(faces)} faces: {dict(cats)}), sliding-elastic two_pass 0, search_radius {SEARCH_RADIUS} mm')


if __name__ == '__main__':
    b135.REFINE_DISTAL = REFINE
    b166.step1()
    with open(os.path.join(RUNS, b166.GEOM, b166.GEOM + '.feb.changes.txt'), 'a', encoding='utf-8') as fo:
        fo.write(f'NOT IN SOURCE: the curves\' last {REFINE[0]} station gaps at the distal end split into {REFINE[1]} along the '
                 f'canal (build_batch135.REFINE_DISTAL; tools/build_batch171.py)\n')
    b166.step2()
    b166.step34()
    p = os.path.join(RUNS, b166.ENDSPR, b166.ENDSPR + '.feb')
    t, msg5 = self_all(open(p, encoding='utf-8').read())
    t, msg6 = wall_damping(t)
    print(msg5)
    write(NAME, t, f'base: claude_diag/runs/{b166.ENDSPR}/{b166.ENDSPR}.feb (tools/build_batch171.py)\n{msg5}\n{msg6}\n')
    print(' '.join(build(171, 'wrapref_endspr_selfall_walls035', {'walls': 0.35}, ['P1', 'P2'], NAME)))
