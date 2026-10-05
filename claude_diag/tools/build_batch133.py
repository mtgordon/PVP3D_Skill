"""Batch 133 (2026-10-01, the user: "make a new version of the non-loft version where the vaginawall is narrowed so that it
does not stick through the levator surface. When the 3mm curved pieces were added it made it go through it from the
start ... make the top half of the curve the avw material and the bottom half the pvw material").
The springs version of the base, L128_tube_r30s_pm (= L91_newline_vwyeoh_rhoi0_seamspr001 + the 3 mm rounded wrap +
SlidingElastic1 seg_up 2 + the PM as a structure), rebuilt with two changes to the wrap (build_batch122's hooks):
  * narrowed: at each wrap station the U's lateral reach is scaled by s (1 = the round 3 mm U; s < 1 squashes it toward
    the walls, the ends on the walls' side faces unchanged), the largest s (steps of 0.05, down to 0.2) for which every
    free wrap node there is at least CLEAR mm in front of the levator's mid-surface (the LA shells are 4 mm thick, so
    2 mm is their surface; CLEAR = 2.5 mm), or more than 8 mm from it; then the minimum over +-2 neighbouring stations, so
    the narrowing is smooth along the canal. In L128_tube_r30s_pm 117 wrap nodes start within 2 mm of the LA mid-surface,
    57 behind it. The walls themselves (AVW, PVW, cervix) are all clear of it already and are not moved;
  * the U's AVW half (segments j < 3 of 6) is Vagina_AVW (domain canal_wrap_hex_avw), its PVW half Vagina_PVW
    (canal_wrap_hex_pvw). Both are the same faithful Yeoh law today, so this changes nothing yet.
  L133_tube_r30s_pm_narrow   built only (the user looks at it first)
usage: py -3.10 build_batch133.py [--check]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS, emit  # noqa: E402
import build_batch122 as b122  # noqa: E402
from build_batch113 import Model113, step2, INNER, OUTER_BOTTOM, closest_on_tri  # noqa: E402
from build_batch132 import seg_up2  # noqa: E402
from febmodel import Feb  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
NAME = 'L133_tube_r30s_pm_narrow'
LA = ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri')
CLEAR, FAR = 2.5, 8.0

f0 = Feb(os.path.join(RUNS, BASE, BASE + '.feb'))
TRI = []
for dm in LA:
    for c in f0.elem_blocks[dm][1].values():
        TRI.append(c[:3])
        if len(c) == 4:
            TRI.append([c[0], c[2], c[3]])
TA = np.array([f0.nodes[t[0]] for t in TRI]); TB = np.array([f0.nodes[t[1]] for t in TRI])
TC = np.array([f0.nodes[t[2]] for t in TRI])
TN = np.cross(TB - TA, TC - TA); TN /= np.linalg.norm(TN, axis=1)[:, None]
TCEN = (TA + TB + TC) / 3


def signed_la(p):
    """(distance to the LA mid-surface, signed + in front of it) at its closest facet among the nearest 12"""
    best = None
    for i in np.argsort(np.linalg.norm(TCEN - p, axis=1))[:12]:
        q, _ = closest_on_tri(p, TA[i], TB[i], TC[i])
        dd = np.linalg.norm(p - q)
        if best is None or dd < best[0]:
            best = (dd, i, q)
    dd, i, q = best
    return dd, (dd if (p - q) @ TN[i] >= 0 else -dd)


def clear(pts):
    for p in pts:
        dd, sd = signed_la(p)
        if dd < FAR and sd < CLEAR:
            return False
    return True


def search(sgn, a0, pts_of_s):
    for s in np.arange(1.0, 0.199, -0.05):
        if clear(pts_of_s(s)):
            return round(float(s), 2)
    return 0.2


def build(model):
    b122.RIN = 3.0
    b122.SPLIT_MAT = True
    # pass 1 (on a scratch model): the raw scale per station
    b122.SCALE_FN, b122.SCALES = search, {}
    scratch = Model113(model.src)
    b122.wrap(scratch)
    raw = dict(b122.SCALES)
    smooth = {}
    for sgn in (-1, 1):
        keys = [k for k in raw if k[0] == sgn]           # in station order
        for i, k in enumerate(keys):
            smooth[k] = min(raw[q] for q in keys[max(0, i - 2):i + 3])
    b122.SCALE_FN, b122.SCALES = (lambda sgn, a0, pts: smooth[(sgn, a0)]), {}
    _, jh, _ = b122.wrap(model)
    assert (jh > 0).all(), 'a wrap hex has a non-positive Jacobian'
    b122.SCALE_FN, b122.SPLIT_MAT = None, False
    seg_up2(model)
    step2(INNER + OUTER_BOTTOM)(model)
    vals = np.array(list(smooth.values()))
    model.log.append(f'NOT IN SOURCE (the wrap narrowed clear of the levator): lateral scale per station min / median / '
                     f'max {vals.min():.2f} / {np.median(vals):.2f} / {vals.max():.2f} ({int((vals < 1).sum())} of '
                     f'{len(vals)} stations narrowed; clearance {CLEAR} mm in front of the LA mid-surface within {FAR} '
                     f'mm, smoothed over +-2 stations); the U\'s AVW half Vagina_AVW (canal_wrap_hex_avw), its PVW half '
                     f'Vagina_PVW (canal_wrap_hex_pvw); hex Jacobians all positive (min {jh.min():.3g})')
    return smooth, jh


if __name__ == '__main__':
    assert not os.path.exists(os.path.join(RUNS, NAME)), f'{NAME} exists; not overwriting'
    m = Model113(os.path.join(RUNS, BASE, BASE + '.feb'))
    smooth, jh = build(m)
    print(m.log[-1])
    for sgn in (-1, 1):
        print(sgn, [smooth[k] for k in smooth if k[0] == sgn])
    if '--check' not in sys.argv:
        emit(NAME, m)
