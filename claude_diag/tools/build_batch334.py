"""Batch 334 (2026-10-19; the user: add the curve-tip fan loft to L288c, the springs line, which has no lofts). Step 1: only the
curve-tip -> PM fan (PM_avw_bottom_left/right_conn + curve_end_PM_conn, 5 single springs a side) as chains woven into one sheet
a side, exactly build_batch262 --fan's recipe for these two families (curve_end_PM_conn first gets its own copy of PM_conn_mat);
every other spring stays single. In L288c the curve ends of these springs are spring-patch anchors (build_batch288), so the
sheet's curve edge lands on the anchors and pulls through the patches. Step 2 (the loft only, fitted): build_batch332.py.
NOT IN SOURCE (chains + sheet).
usage: py -3.10 build_batch334.py SRC_RUN NEW_RUN"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch111 import Model111, SHEETS  # noqa: E402
from build_batch262 import own_curve_end_mat  # noqa: E402
from variants2 import renumber_discrete_materials  # noqa: E402
from variants import emit  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def main(src, new):
    m = Model111(os.path.join(RUNS, src, src + '.feb'))
    m.log.append(f'{src} -> the curve-tip fan only as chains in a sheet (build_batch334.py)')
    own_curve_end_mat(m)
    SHEETS['PM_avw_bottom_left'] = ('PM_avw_bottom_left_fan', ['PM_avw_bottom_left_conn', 'curve_end_PM_conn'], 'x-')
    SHEETS['PM_avw_bottom_right'] = ('PM_avw_bottom_right_fan', ['PM_avw_bottom_right_conn', 'curve_end_PM_conn'], 'x+')
    gone = set()
    for sheet in ('PM_avw_bottom_left', 'PM_avw_bottom_right'):
        gone.update(m.chains(sheet, 0.05, 1.0))
    disc = m.root.find('Discrete')
    for s in sorted(gone):
        pairs, mat, b, ds = m._spring_set(s)
        disc.remove(b); disc.remove(mat); m.mesh.remove(ds)
        m.log.append(f'{s}: its {len(pairs)} single springs removed (now chains)')
    renumber_discrete_materials(m)
    assert not os.path.exists(os.path.join(RUNS, new)), new
    emit(new, m)
    print('\n'.join(m.log))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
