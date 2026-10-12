"""Batch 162 (2026-10-05, the user: connect the PM-AVW-bottom connectors (E28102 and the others, both connectors on both sides)
"to the new curved edge of the tube instead of the original AVW part. Ideally one of the connectors on each side would connect
to the AVW and the other to the pvw").
The newest tube geometry L143_thin275_narrow3_wrapend_noflare has its 27 PM_conn-family springs on the PM's inner arc (built
before the outer-arc rule of 10-02), so it is first rebuilt with them on the outer arc, as L149: L143's own stage file before
its PM_conn step (scratch_build/L143_..._tmppre) + build_batch149.outer_arc + the curves' inner faces in the lumen contact
(SlidingElastic1Primary / Secondary copied from L143, the last step of build_batch143).
The 4 PM_avw_bottom springs (source CONN3D2 246-249, law PM-AVW_bottom_LA-x%stiff = PM-LA x 1e-5) run from their anchors on
the PM (tied to the PM shell) to the AVW's side edge, 11.1-12.2 mm (8.2 / 9.0 mm in the source: the 3 mm narrowing moved the
AVW ends inward). Each is reattached to the curve (canal_wrap_hex_*) on its side, at the curve's surface node nearest the
spring's PM end (the user's choice): the upper spring of each side (left_1 / right_1) to the AVW half, the lower (left_2 /
right_2) to the PVW half; nodes in a wall domain, in both halves, or in a linear constraint are excluded. Law unchanged.
  L162_thin275_wrapend_outer                 L143 with PM_conn on the outer arc (the reference)
  L162_thin275_wrapend_outer_avwbcurve       + the 4 springs moved to the curves, source law (near zero)
  L162_thin275_wrapend_outer_avwbcurve_full  + the 4 springs moved, with the source's PM-LA-x%stiff law at full (x ~1e5)
usage: py -3.10 build_batch162.py"""
import os
import re
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch149 import outer_arc  # noqa: E402
from febmodel import Feb  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SCRATCH = os.path.join(JOBS, 'claude_diag', 'scratch_build')
L143 = 'L143_thin275_narrow3_wrapend_noflare'
PRE = os.path.join(SCRATCH, L143 + '_tmppre', L143 + '_tmppre.feb')
SRC143 = os.path.join(RUNS, L143, L143 + '.feb')
INP = os.path.join(JOBS, 'PVP3DModel_job.inp')
REF = 'L162_thin275_wrapend_outer'
SETS = ('PM_avw_bottom_left_conn', 'PM_avw_bottom_right_conn')
MATS = ('PM_avw_bottom_left_conn_mat', 'PM_avw_bottom_right_conn_mat')
WALLS = ('_PickedSet347', '_PickedSet346', '_PickedSet64', '_PickedSet66')
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def block(t, tag, name):
    m = list(re.finditer(r'<%s [^>]*name="%s"' % (tag, re.escape(name)), t))
    assert len(m) == 1, (tag, name, len(m))
    i = m[0].start()
    return i, t.index(f'</{tag}>', i) + len(f'</{tag}>')


def copy_surface(t, src, name):
    i, j = block(t, 'Surface', name)
    a, b = block(src, 'Surface', name)
    return t[:i] + src[a:b] + t[j:]


def elems(t, name):
    i, j = block(t, 'Elements', name)
    return [[int(v) for v in c.split(',')] for c in re.findall(r'<elem id="\d+">([^<]+)</elem>', t[i:j])]


def move_springs(t, path):
    f = Feb(path)
    X = {n: np.array(p, float) for n, p in f.nodes.items()}
    bnd, mem = {}, {}
    for name, tag in (('canal_wrap_hex_avw', 'AVW'), ('canal_wrap_hex_pvw', 'PVW')):
        E = elems(t, name)
        cnt = Counter(tuple(sorted(e[i] for i in fc)) for e in E for fc in HF)
        bnd[tag] = {n for fc, c in cnt.items() if c == 1 for n in fc}
        mem[tag] = {n for e in E for n in e}
    wall = {n for w in WALLS for e in elems(t, w) for n in e}
    constrained = {int(v) for v in re.findall(r'<node id="(\d+)" bc="[xyz]">', t)}
    cand = {'AVW': bnd['AVW'] - mem['PVW'] - wall - constrained, 'PVW': bnd['PVW'] - mem['AVW'] - wall - constrained}
    rep = []
    for nm in SETS:
        i, j = block(t, 'DiscreteSet', nm)
        blk = t[i:j]
        pairs = [tuple(map(int, d.split(','))) for d in re.findall(r'<delem>([^<]+)</delem>', blk)]
        assert len(pairs) == 2, (nm, pairs)
        for (anc, old), half in zip(pairs, ('AVW', 'PVW')):        # source order: _1 (upper), _2 (lower)
            new = min(cand[half], key=lambda n: np.linalg.norm(X[n] - X[anc]))
            blk = blk.replace(f'<delem>{anc},{old}</delem>', f'<delem>{anc},{new}</delem>')
            rep.append(f'{nm} {anc}->{old} ({np.linalg.norm(X[anc] - X[old]):.2f} mm) => {anc}->{new} on the curve\'s {half} '
                       f'half ({np.linalg.norm(X[anc] - X[new]):.2f} mm; node at {", ".join(f"{v:.2f}" for v in X[new])})')
        t = t[:i] + blk + t[j:]
    return t, ('NOT IN SOURCE: the 4 PM_avw_bottom springs reattached from the AVW\'s side edge to the curves (canal_wrap_hex_*), '
               'at the curve surface node nearest each spring\'s PM end, upper spring to the AVW half, lower to the PVW half; '
               'law unchanged: ' + '; '.join(rep))


def full_law(t):
    s = open(INP, encoding='ISO-8859-1').read()
    i = s.index('*Connector Behavior, name=PM-LA-x%stiff\n')
    rows = re.findall(r'^\s*([-\d.eE+]+),\s*([-\d.eE+]+)\s*$', s[i:s.index('*Connector Behavior', i + 10)], re.M)
    pts = ''.join(f'\n\t\t\t\t\t<pt>{float(u):.6g},{float(F):.6g}</pt>' for F, u in rows)    # Abaqus: force first
    assert len(rows) == 11
    for mat in MATS:
        i = t.index(f'name="{mat}" type="nonlinear spring">')
        a = t.index('<points>', i) + len('<points>')
        b = t.index('</points>', a)
        t = t[:a] + pts + '\n\t\t\t\t' + t[b:]
    return t, (f'NOT IN SOURCE: {", ".join(MATS)} -> the source\'s PM-LA-x%stiff law at full (elongation as in the source, '
               f'{rows[1][0]} N at {rows[1][1]} mm ... {rows[-1][0]} N at {rows[-1][1]} mm; the source\'s PM-AVW_bottom law is it '
               f'x ~1e-5)')


def write(name, txt, notes):
    d = os.path.join(RUNS, name)
    assert not os.path.exists(d), f'{name} exists; not overwriting'
    os.makedirs(d)
    p = os.path.join(d, name + '.feb')
    open(p, 'w', encoding='utf-8', newline='').write(txt)
    with open(p + '.changes.txt', 'w', encoding='utf-8') as fo:
        fo.write(notes)
    print('wrote', p)
    return p


if __name__ == '__main__':
    t0 = open(PRE, encoding='utf-8').read()
    t143 = open(SRC143, encoding='utf-8').read()
    pre_log = open(PRE + '.changes.txt', encoding='utf-8').read()
    t, msg = outer_arc(t0, PRE)
    for s in ('SlidingElastic1Primary', 'SlidingElastic1Secondary'):
        t = copy_surface(t, t143, s)
    note = (f'{REF} (tools/build_batch162.py) from {os.path.relpath(PRE, JOBS)} (= build_batch143 before its PM_conn step):\n'
            f'{pre_log}{msg}\nthe curves\' inner faces in the lumen contact: SlidingElastic1Primary / Secondary as in {L143}\n')
    print(msg)
    p_ref = write(REF, t, note)
    t1, msg1 = move_springs(t, p_ref)
    print(msg1)
    write(REF + '_avwbcurve', t1, f'base: claude_diag/runs/{REF}/{REF}.feb (tools/build_batch162.py)\n{msg1}\n')
    t2, msg2 = full_law(t1)
    print(msg2)
    write(REF + '_avwbcurve_full', t2, f'base: claude_diag/runs/{REF}_avwbcurve/{REF}_avwbcurve.feb (tools/build_batch162.py)\n'
                                       f'{msg1}\n{msg2}\n')
