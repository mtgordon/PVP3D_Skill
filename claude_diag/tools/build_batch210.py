"""Batch 210 (2026-10-07 night, unattended; the user: "what is prompting the large difference between this model and the
Abaqus model?"). The audit: on the source-support lips model L209_src_sl01_P1 (source support values, contact-only seam with
x0.1 lips), undo one of the model's other departures from the source per variant:
  chainmass   the arcus chain mass at the Abaqus truss density (chain_mass_mat 0.0011 -> 0.00011; the line has x10)
  nodamp      no added mass damping on the walls / lips / PM (wall_damping_*, pm_damping removed; settle_damping acts only
              after t 1)
  loadla      Load-LA at the source's 0.014 MPa (the line has 1/3, the user's earlier choice)
usage: py -3.10 build_batch210.py SRC_RUN MODE NEW_NAME"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')

if __name__ == '__main__':
    src, mode, name = sys.argv[1:4]
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    if mode == 'chainmass':
        t, n = re.subn(r'(<material id="\d+" name="chain_mass_mat" type="linear truss">\s*<density>)0\.0011(</density>)',
                       r'\g<1>0.00011\2', t)
        assert n == 1
        msg = 'as the source: the arcus chain mass at the Abaqus truss density (chain_mass_mat density 0.0011 -> 0.00011)'
    elif mode == 'nodamp':
        names = re.findall(r'<body_load name="((?:wall_damping_[^"]+)|pm_damping)" type="mass damping"', t)
        for nm in names:
            t, n = re.subn(r'\t*<body_load name="%s" type="mass damping"[^>]*>.*?</body_load>\n' % re.escape(nm), '', t,
                           flags=re.S)
            assert n == 1, nm
        msg = f'as the source: the added mass damping removed ({", ".join(names)}); settle_damping (after t 1) kept'
    elif mode == 'loadla':
        t, n = re.subn(r'(<surface_load name="Load-LA" surface="Load-LA_surf" type="pressure">\s*<pressure lc="1">)'
                       r'[^<]+(</pressure>)', r'\g<1>0.014\2', t)
        assert n == 1
        msg = 'as the source: Load-LA 0.014 MPa (was 1/3 of it, the user\'s earlier choice)'
    else:
        raise SystemExit('mode: chainmass | nodamp | loadla')
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch210.py)\n{msg}\n')
    print(msg)
