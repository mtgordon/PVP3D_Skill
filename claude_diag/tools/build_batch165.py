"""Batch 165 (2026-10-05, the user: the runs "tend to be waving once they get loaded. Can we add more mass or dampening in to
reduce that?" ... "do it to the 35% version as that is where I'm setting it more").
Measured on L164_avwbend_fixed_self_walls07_P2 (t 0.95-1): the walls oscillate at ~30 Hz (wall mean velocity +35 / -34 mm/s,
period ~0.035 s); rms about the trend: cervix 0.44 mm, AVW 0.15, PVW 0.13, curves 0.11-0.12, PeB 0.03, LA <= 0.016, chain
0.003. The existing mass damping (settle_damping, C 20/s, every mass) only switches on after t = 1, where the runs end.
NOT IN SOURCE: mass-proportional damping on the walls only (one body_load "mass damping" per domain, elem_set = the domain;
verified on a two-block test that elem_set confines it), C = 200/s from t = 0 (no load curve): ~0.55 of critical for the
30 Hz mode, while its lag on the wall's slow response (C / omega^2 ~ 6 ms) and its force (C rho v ~ 0.1 % of the pressure
load) stay small. Domains: _PickedSet347 (AVW), _PickedSet346 (cervix), _PickedSet64 (PVW), canal_wrap_hex_avw / _pvw.
  L165_self_walls035_wdamp200_P1   L164_avwbend_fixed_self_walls035_P1 + the wall damping
  L165_self_walls035_wdamp200_P2   L164_avwbend_fixed_self_walls035_P2 + the wall damping
usage: py -3.10 build_batch165.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
C = 200
DOMAINS = ('_PickedSet347', '_PickedSet346', '_PickedSet64', 'canal_wrap_hex_avw', 'canal_wrap_hex_pvw')
BUILDS = (('L165_self_walls035_wdamp200_P1', 'L164_avwbend_fixed_self_walls035_P1'),
          ('L165_self_walls035_wdamp200_P2', 'L164_avwbend_fixed_self_walls035_P2'))


def wall_damping(t, domains=DOMAINS):
    for d in domains:
        assert t.count(f'<SolidDomain name="{d}"') == 1, d
    add = ''.join(f'\t\t<body_load name="wall_damping_{d.strip("_")}" type="mass damping" elem_set="{d}">\n'
                  f'\t\t\t<C>{C}</C>\n\t\t</body_load>\n' for d in domains)
    i = t.index('\t</Loads>')
    return t[:i] + add + t[i:], (f'NOT IN SOURCE: mass damping on the walls only, C = {C}/s from t = 0 (body_load "mass damping", '
                                 f'elem_set {", ".join(domains)}), to damp the walls\' ~30 Hz waving once loaded')


if __name__ == '__main__':
    for name, src in BUILDS:
        p = os.path.join(RUNS, src, src + '.feb')
        t, msg = wall_damping(open(p, encoding='utf-8').read())
        write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch165.py)\n{msg}\n')
