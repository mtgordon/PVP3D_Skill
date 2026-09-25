"""Convert an Abaqus *Hyperelastic, marlow + *Uniaxial Test Data material to FEBio's uncoupled Yeoh.

Why Yeoh, not a uniaxial Ogden fit: Marlow's deviatoric energy depends on I1 only, so its biaxial and planar
response follow from the uniaxial data. Yeoh is also I1-only, so a good uniaxial fit reproduces Marlow in
every stretch state. A 1-term Ogden fitted to the same uniaxial data matched to +-10 % in uniaxial tension
but was 2-3x softer in equibiaxial stretch at 20-40 % strain (febio-xml-format.md gotcha 23).

usage: py -3.10 fit_yeoh.py data.txt [--nmax 4] [--nu 0.47] [--ogden c1,m1[,c2,m2,...]]
  data.txt : one "stress strain" pair per line (comma or whitespace), Abaqus *Uniaxial Test Data order
             (nominal stress FIRST, then nominal strain). Check the column order against the source.
  --ogden  : FEBio Ogden c_i,m_i to compare against (e.g. an earlier fit), incompressible.
Output: Yeoh c_i per N (relative-error least squares, linear in c_i: unique, no multi-start), fit errors,
stability (all c_i >= 0 guarantees W'>0, W''>=0), the FEBio <material> block, and uniaxial / planar /
equibiaxial nominal stress vs the Marlow construction. Also k: from mu0 (matches nu at zero strain only) and
the k that reproduces Abaqus' constant-nu volume change J = lambda^(1-2nu) at 20-50 % uniaxial strain
(FEBio's constant k otherwise lets the material lose volume stiffness as it stiffens; gotcha 23).
"""
import argparse
import re

import numpy as np


def load(path):
    rows = [list(map(float, re.split(r'[,\s]+', ln.strip())))[:2] for ln in open(path) if re.match(r'\s*[-\d.]', ln)]
    s, e = np.array(rows).T
    order = np.argsort(e)
    return e[order], s[order]


def i1_states():
    return {'uniaxial': (lambda l: (l, l ** -0.5, l ** -0.5), lambda l, w: 2 * w * (l - l ** -2)),
            'planar': (lambda l: (l, 1.0, 1.0 / l), lambda l, w: 2 * w * (l - l ** -3)),
            'equibiaxial': (lambda l: (l, l, l ** -2), lambda l, w: 2 * w * (l - l ** -5))}


def marlow_dw(I1, eps, sig):
    lo, hi = 1.0, 10.0          # uniaxial stretch with the same I1 (monotonic for l >= 1)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mid * mid + 2 / mid < I1 else (lo, mid)
    lu = 0.5 * (lo + hi)
    return np.interp(lu - 1, eps, sig) / (2 * (lu - lu ** -2))


def ogden_T(state, l, cm):
    out = 0.0
    for c, m in cm:
        tail = {'uniaxial': -m / 2 - 1, 'planar': -m - 1, 'equibiaxial': -2 * m - 1}[state]
        out += c / m * (l ** (m - 1) - l ** tail)
    return out


def fit(eps, sig, n):
    keep = eps > 0
    l = 1 + eps[keep]
    I1 = l ** 2 + 2 / l
    A = np.array([2 * (l - l ** -2) * i * (I1 - 3) ** (i - 1) for i in range(1, n + 1)]).T
    c, *_ = np.linalg.lstsq(A / sig[keep, None], np.ones(keep.sum()), rcond=None)
    return c, (A @ c - sig[keep]) / sig[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data')
    ap.add_argument('--nmax', type=int, default=4)
    ap.add_argument('--nu', type=float, default=0.47)
    ap.add_argument('--ogden', default='')
    a = ap.parse_args()
    eps, sig = load(a.data)
    best = None
    for n in range(1, a.nmax + 1):
        c, err = fit(eps, sig, n)
        stable = bool(np.all(c >= 0))
        print(f'N={n}: c = {", ".join("%.6g" % v for v in c)} | max |rel err| {np.abs(err).max():.1%} | '
              f'all c_i >= 0 (stable): {stable}')
        if best is None and stable and np.abs(err).max() <= 0.10:
            best = (n, c)
    n, c = best or (2, fit(eps, sig, 2)[0])
    mu0 = 2 * c[0]
    k0 = 2 * mu0 * (1 + a.nu) / (3 * (1 - 2 * a.nu))
    print(f'\nchosen N={n} (smallest stable fit within 10 %); mu0 = 2 c1 = {mu0:.5g}; k from mu0, nu={a.nu}: {k0:.4g}')
    ks = []
    for e in (0.2, 0.3, 0.4, 0.5):
        l = 1 + e
        J = l ** (1 - 2 * a.nu)
        ks.append((e, np.interp(e, eps, sig) * l / 3 * J / np.log(J)))
    print('k reproducing Abaqus constant-nu J = lambda^(1-2nu) in uniaxial tension: '
          + ', '.join(f'{k:.3g} at {e:.0%}' for e, k in ks))
    print('\nFEBio block (choose k for the working strain range):')
    print('<material id="?" name="?" type="Yeoh">\n  <density>?</density>\n'
          + ''.join(f'  <c{i}>{v:.7g}</c{i}>\n' for i, v in enumerate(c, 1)) + f'  <k>{k0:.4g}</k>\n</material>')
    cm = list(zip(*[iter(map(float, a.ogden.split(',')))] * 2)) if a.ogden else []
    dw = lambda I1: sum((i + 1) * ci * (I1 - 3) ** i for i, ci in enumerate(c))
    for state, (st, t_i1) in i1_states().items():
        print(f'\n{state}: nominal stress (incompressible)  stretch | Marlow | Yeoh' + (' | Ogden' if cm else ''))
        for l in (1.1, 1.2, 1.3, 1.4, 1.5):
            I1 = sum(v ** 2 for v in st(l))
            row = f'  {l:4.2f} | {t_i1(l, marlow_dw(I1, eps, sig)):.5g} | {t_i1(l, dw(I1)):.5g}'
            print(row + (f' | {ogden_T(state, l, cm):.5g}' if cm else ''))


if __name__ == '__main__':
    main()
