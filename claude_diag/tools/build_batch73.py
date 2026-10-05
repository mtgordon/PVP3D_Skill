"""Batch 73 (2026-09-27, Test A of NEXT_SESSION_PROMPT_2026-09-26b.md): which fitted loft makes the lofts line faster.

The lofts line (L26_lofts_la3_pm: 308 steps, 3068 Newton iterations, 21:20) beats the springs line (L26_springs_la3_pm:
419 steps, 4693 iterations, 23:14) with a third fewer iterations, although its 50050 equations (the loft shells) against
33601 make each iteration ~40 % dearer. One change each from L26_springs_la3_pm: one of the user's four fitted loft
families put back in place of its Abaqus connectors (variants7.Model7.connectors_to_lofts, copied from L26_lofts_la3_pm;
all four back reproduces L26_lofts_la3_pm exactly), plus fresh copies of both bases as controls, to run in one batch at one
thread count.
  L66_springs_la3_pm_avwparaloft  AVW-Para-L/R lofts (+ AVW-Para-L's two tied-node-on-facet ties); 35 + 35 springs out
  L66_springs_la3_pm_clloft       CL-L/R lofts; CL-L/R_conn_1..3 out; BC-CL-Left/Right back to 65 nodes
  L66_springs_la3_pm_uslloft      USL-L/R lofts (+ USL-L's two tied-node-on-facet ties); USL-L/R_conn_1..2 out;
                                  BC-USL-Left/Right back to 41 nodes
  L66_springs_la3_pm_pmloft       the PM loft; PM_conn out
  L66_springs_la3_pm_ctrl         identical copy of L26_springs_la3_pm
  L66_lofts_la3_pm_ctrl           identical copy of L26_lofts_la3_pm
usage: py -3.10 build_batch73.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants7 import Model7, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SPRINGS = os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')
LOFTS = os.path.join(RUNS, 'L26_lofts_la3_pm', 'L26_lofts_la3_pm.feb')
WANT = set(sys.argv[1:])


def loft_back(family):
    def fn(m):
        m.connectors_to_lofts(LOFTS, family)
    return fn


def control(m):
    m.log.append('identical copy (the control for batch 73: run in the same batch at the same thread count as the '
                 'one-loft variants)')


BUILDS = (('L66_springs_la3_pm_avwparaloft', SPRINGS, loft_back('AVW-Para')),
          ('L66_springs_la3_pm_clloft', SPRINGS, loft_back('CL')),
          ('L66_springs_la3_pm_uslloft', SPRINGS, loft_back('USL')),
          ('L66_springs_la3_pm_pmloft', SPRINGS, loft_back('PM')),
          ('L66_springs_la3_pm_ctrl', SPRINGS, control),
          ('L66_lofts_la3_pm_ctrl', LOFTS, control))
if __name__ == '__main__':
    for name, base, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
        m = Model7(base)
        fn(m)
        emit(name, m)
