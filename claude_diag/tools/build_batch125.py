"""Batches 125-126 (2026-09-29 ~22:25, the user: "Look into why the tube's P1 case was slow, before adopting it. Then try
adjusting the vaginal wall stiffness to try and get closer to the paper numbers").

Batch 125: the slow tube P1 (L124_tube_round_r30_paperP1: t = 1 but 114 failed steps, 4:22, all at t 0.965-1.0). From
t ~0.965 the AVW / cervix lumen-edge seam nodes 2076-2079 / 2152-2155 (both sides, y -9 .. -14) chatter on the PVW's edge
facets: with the rounded end they hang 3 mm below the wrap's attachment, held only by the contact SlidingElastic1
(seg_up 0 = unlimited segment updates). The failed steps end with the line search at 1e-6 ("Zero linestep size") while
the energy and displacement norms are met. One change each, solver / contact settings only (should not move the answer):
  L125_tube_round_r30_paperP1_lstol0     line search off (lstol 0.9 -> 0)
  L125_tube_round_r30_paperP1_segup2     SlidingElastic1 seg_up 0 -> 2 (the PVW_LA contact's value)
  L125_tube_round_r30_paperP1_lsjac0     ls_check_jacobians 1 -> 0 (the line search no longer shrinks a step that
                                         inverts an element in a trial state)
Batch 126: the vaginal walls' stiffness on the tube, both paper cases (Vagina_AVW, Vagina_PVW, Vagina_Cervix Yeoh c1, c2, k
x F; the wrap uses Vagina_AVW, so it follows the anterior factor). NOT FAITHFUL (the faithful Yeoh refit is x1).
  L126_tube_r30_walls{tag}_paperP1 / P2
usage: py -3.10 build_batch125.py NAME [NAME ...]   (see VARIANTS; an existing run folder is never overwritten)"""
import os
import re
import sys

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
C1, C2, K = 0.005964127, 0.01883547, 1.0


def sub1(t, pat, rep, label):
    t2, n = re.subn(pat, rep, t, count=1, flags=re.S)
    assert n == 1, label
    return t2


def contact_set(name, tag, old, new):
    def fn(t):
        i = t.index(f'<contact name="{name}"')
        j = t.index('</contact>', i)
        blk = t[i:j]
        a, b = f'<{tag}>{old}</{tag}>', f'<{tag}>{new}</{tag}>'
        assert blk.count(a) == 1, (name, tag)
        return t[:i] + blk.replace(a, b) + t[j:], f'contact {name} {tag} {old} -> {new}'
    return fn


def solver_set(tag, old, new):
    def fn(t):
        return sub1(t, f'<{tag}>{old}</{tag}>', f'<{tag}>{new}</{tag}>', tag), f'solver {tag} {old} -> {new}'
    return fn


def walls(fa, fp):
    """anterior (Vagina_AVW, Vagina_Cervix; the wrap is Vagina_AVW) x fa, posterior (Vagina_PVW) x fp"""
    def fn(t):
        for nm, f in (('Vagina_AVW', fa), ('Vagina_Cervix', fa), ('Vagina_PVW', fp)):
            i = t.index(f'name="{nm}" type="Yeoh">')
            j = t.index('</material>', i)
            blk = t[i:j]
            for tag, v in (('c1', C1), ('c2', C2), ('k', K)):
                blk = sub1(blk, f'<{tag}>[^<]*</{tag}>', f'<{tag}>{v * f:.7g}</{tag}>', nm + tag)
            t = t[:i] + blk + t[j:]
        return t, (f'NOT FAITHFUL: vaginal walls Yeoh c1, c2, k x {fa} (Vagina_AVW, Vagina_Cervix; the wrap is Vagina_AVW) '
                   f'and x {fp} (Vagina_PVW)')
    return fn


T1, T2 = 'L124_tube_round_r30_paperP1', 'L124_tube_round_r30_paperP2'
VARIANTS = {
    'L125_tube_round_r30_paperP1_lstol0': (T1, solver_set('lstol', '0.9', '0')),
    'L125_tube_round_r30_paperP1_segup2': (T1, contact_set('SlidingElastic1', 'seg_up', '0', '2')),
    'L125_tube_round_r30_paperP1_lsjac0': (T1, solver_set('ls_check_jacobians', '1', '0')),
}
for tag, fa, fp in (('x2', 2, 2), ('x05', 0.5, 0.5), ('a2p05', 2, 0.5), ('a2', 2, 1), ('p05', 1, 0.5),
                    ('x15', 1.5, 1.5), ('x3', 3, 3), ('a3p05', 3, 0.5), ('a2p07', 2, 0.7), ('a4', 4, 1), ('a4p05', 4, 0.5)):
    for c, src in (('P1', T1), ('P2', T2)):
        VARIANTS[f'L126_tube_r30_walls{tag}_paper{c}'] = (src, walls(fa, fp))

# Batch 127: as 126, on the tube with the P1 speed fix (SlidingElastic1 seg_up 2; L125_tube_round_r30_paperP1_segup2 gave the
# same answer as L124_tube_round_r30_paperP1 in 0:49 instead of 4:22)
def both(f1, f2):
    def fn(t):
        t, a = f1(t)
        t, b = f2(t)
        return t, a + '; ' + b
    return fn


SEG = contact_set('SlidingElastic1', 'seg_up', '0', '2')
for tag, fa, fp in (('a2p05', 2, 0.5), ('a3p05', 3, 0.5), ('a3p03', 3, 0.3), ('a4p05', 4, 0.5), ('a4p03', 4, 0.3),
                    ('a2p03', 2, 0.3), ('a3', 3, 1), ('a4', 4, 1), ('p03', 1, 0.3), ('a15p05', 1.5, 0.5),
                    ('a3p07', 3, 0.7), ('a2p07', 2, 0.7), ('x2', 2, 2), ('a5p05', 5, 0.5), ('a6p05', 6, 0.5), ('x05', 0.5, 0.5), ('a2', 2, 1), ('p05', 1, 0.5), ('p02', 1, 0.2), ('a2p02', 2, 0.2)):
    for c, src in (('P1', T1), ('P2', T2)):
        VARIANTS[f'L127_tube_r30s_walls{tag}_paper{c}'] = (src, both(walls(fa, fp), SEG))

# the healthy (non-paper) tube with the same walls, to show how far each wall change moves the base case
for tag, fa, fp in (('a2', 2, 1), ('a2p05', 2, 0.5), ('a2p03', 2, 0.3), ('p03', 1, 0.3)):
    VARIANTS[f'L127_tube_r30s_walls{tag}'] = ('L123_tube_round_r30', both(walls(fa, fp), SEG))
VARIANTS['L127_tube_r30s'] = ('L123_tube_round_r30', SEG)

if __name__ == '__main__':
    for name in sys.argv[1:]:
        src, fn = VARIANTS[name]
        d = os.path.join(RUNS, name)
        assert not os.path.exists(d), f'{name} exists; not overwriting'
        t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
        t2, what = fn(t)
        assert t2 != t
        os.makedirs(d)
        open(os.path.join(d, name + '.feb'), 'w', encoding='utf-8', newline='').write(t2)
        open(os.path.join(d, name + '.feb.changes.txt'), 'w', encoding='utf-8').write(f'{name}: {src} + {what}\n')
        print(name, ':', what)
