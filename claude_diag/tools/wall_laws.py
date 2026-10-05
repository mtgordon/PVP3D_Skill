"""The vaginal-wall laws side by side (2026-09-28, batch 99): the source's Marlow data, the faithful Yeoh refit, the Yeoh
refit with c1 scaled (x1.25, x1.5, x2, x3: NOT FAITHFUL) and the converter's Ogden walls, in uniaxial, planar and equibiaxial
stretch (incompressible nominal stress; the skill's scripts/fit_yeoh.py machinery).
usage: py -3.10 wall_laws.py [--ratio]   (--ratio: each law relative to the source's Marlow law instead of kPa)
"""
import os
import sys

sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from fit_yeoh import load, i1_states, marlow_dw, ogden_T  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scratch_2026-09-27', 'vagina_marlow.txt')
C1, C2 = 0.005964127, 0.01883547
LAWS = [('Yeoh (faithful)', ('yeoh', C1, C2)), ('c1 x1.25', ('yeoh', 1.25 * C1, C2)), ('c1 x1.5', ('yeoh', 1.5 * C1, C2)),
        ('c1 x2', ('yeoh', 2 * C1, C2)), ('c1 x3', ('yeoh', 3 * C1, C2)), ('Ogden AVW', ('ogden', [(0.1, 1), (0.02, 5.159)])),
        ('Ogden PVW/cx', ('ogden', [(0.0887, 1), (0.0193, 5.159)]))]


def stress(law, state, l, st, t_i1):
    if law[0] == 'yeoh':
        I1 = sum(v ** 2 for v in st(l))
        return t_i1(l, law[1] + 2 * law[2] * (I1 - 3))
    return ogden_T(state, l, law[1])


def main():
    ratio = '--ratio' in sys.argv
    eps, sig = load(DATA)
    print('initial shear modulus [kPa]: ' + ', '.join(
        f'{n} {1e3 * (2 * p[1] if p[0] == "yeoh" else sum(c for c, m in p[1]) / 2):.1f}' for n, p in LAWS))
    for state, (st, t_i1) in i1_states().items():
        print(f'\n{state}: nominal stress ' + ('/ Marlow' if ratio else '[kPa]'))
        print('stretch | Marlow | ' + ' | '.join(n for n, _ in LAWS))
        for l in (1.05, 1.1, 1.2, 1.3, 1.5):
            I1 = sum(v ** 2 for v in st(l))
            mar = t_i1(l, marlow_dw(I1, eps, sig))
            vals = [stress(p, state, l, st, t_i1) for _, p in LAWS]
            if ratio:
                print(f'{l:4.2f} | {1e3 * mar:.2f} kPa | ' + ' | '.join(f'{v / mar:.2f}' for v in vals))
            else:
                print(f'{l:4.2f} | {1e3 * mar:.2f} | ' + ' | '.join(f'{1e3 * v:.2f}' for v in vals))


if __name__ == '__main__':
    main()
