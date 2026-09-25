"""Batch 24 (2026-09-24, user: "I like the lofts and Edge BC on it"; "make the CL/USL tissue properties similar to
Abaqus, accounting for the number of connectors it had, compared with the cross-sectional area it has now").
All with the line search on, to t = 1.0.
Fast base:
  L9_fitall_edge   L8_fitall (every loft fitted, strip model) + the PM_PeB anchor edges fully fixed: the new reference
  L10_sec          L9_fitall_edge + the CL/USL lofts' material from the cross-section (cl_usl_section.py: sum of the
                   connectors' forces / (average width x 0.49 mm) at strain u / mean length, 2-term Ogden, rms 9 %;
                   k = 250 mu0), density kept (7.8e-07, as L8_fitall)
  L10_secrho       L10_sec + the CL/USL lofts at tissue density (1.06e-9; at 7.8e-07 each CL loft weighs ~250 g and each
                   USL loft ~150 g; all the tissue weighs 69 g; the connectors they replace are massless)
Like-for-like base (Abaqus chain mass):
  L9D_fitall_edge  L8D_fitall + the PM_PeB anchor edges fully fixed
  L10D_sec         L9D_fitall_edge + the CL/USL cross-section material
usage: py -3.10 build_batch24.py [NAME ...]   (default: all; an existing run folder is never overwritten)
"""
import os
import sys

from variants6 import Model6, JOBS, emit
from cl_usl_section import section_fits

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
TISSUE = 1.06e-9
CLUSL = ('CL-L_fan', 'CL-R_fan', 'USL-L_fan', 'USL-R_fan')
WANT = set(sys.argv[1:])
SEC = {}


def want(name):
    return not WANT or name in WANT


def new(name, base):
    assert not os.path.exists(os.path.join(RUNS, name)), f'{name} exists; not overwriting'
    return Model6(os.path.join(RUNS, base, base + '.feb'))


def edge(m):
    for loft in ('PM_PeB_Left_fan', 'PM_PeB_Right_fan'):
        m.fix_anchor_edge(loft, 'BC-RP-PM-PeB-origins')


def section(m):
    if not SEC:
        SEC.update(section_fits(m.src))
    for loft in CLUSL:
        r = SEC[loft]
        m.set_loft_ogden(loft, r['t2'], 250.0, note=(
            f"cross-section fit (cl_usl_section.py): {len(r['rows'])} Abaqus connectors, total force / (average width "
            f"{r['w']:.2f} mm x 0.49 = {r['A']:.2f} mm^2) at strain u / {r['Lm']:.1f} mm, 2-term Ogden rms "
            f"{100 * r['e2']:.0f} %"))


if want('L9_fitall_edge'):
    m = new('L9_fitall_edge', 'L8_fitall'); edge(m); emit('L9_fitall_edge', m)
if want('L10_sec'):
    m = new('L10_sec', 'L8_fitall'); edge(m); section(m); emit('L10_sec', m)
if want('L10_secrho'):
    m = new('L10_secrho', 'L8_fitall'); edge(m); section(m); m.set_lofts_density(CLUSL, TISSUE); emit('L10_secrho', m)
if want('L9D_fitall_edge'):
    m = new('L9D_fitall_edge', 'L8D_fitall'); edge(m); emit('L9D_fitall_edge', m)
if want('L10D_sec'):
    m = new('L10D_sec', 'L8D_fitall'); edge(m); section(m); emit('L10D_sec', m)
