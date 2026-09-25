"""Fit a 1-term Ogden membrane material for the P-arcus fan so that the fan pulls like the 13 connectors it replaces.

Model: connector i (length L_i, from arcus node to vaginal-wall node) becomes a fan strip of width w_i (mean of its
spacing along the vaginal-wall edge and along the arcus edge) and thickness t. Pulling the arcus edge a uniform
distance u away from the wall gives stretch lambda_i = 1 + u/L_i. The fan force is sum_i P(lambda_i) w_i t, with
uniaxial incompressible Ogden P = (c/m)(lambda^(m-1) - lambda^(-m/2-1)) (FEBio c = 2 mu, m = alpha). The target is
sum_i F_conn(u) (Abaqus LA-Y-parcus table, force-first). Least squares over u in [2, 40] mm, relative errors.
usage: py -3.10 fit_fan_to_connectors.py
"""
import numpy as np

from febmodel import Feb
import os
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import RUNS_DIR

Fk = np.array([0, .11667, .28, .56, .84, 1.26, 1.75, 2.42667, 3.26667, 4.41, 5.92667])
uk = np.array([0, 4.7, 9.515, 13.985, 18.8, 23.5, 28.2, 33.015, 37.715, 42.415, 47.115])
PVW = [1847, 1913, 1912, 1911, 1910, 1909, 1908, 1907, 1906, 1905, 1904, 1902, 1900]
t = 0.125
f = Feb(os.path.join(RUNS_DIR, 'LP2_yeoh_k1', 'LP2_yeoh_k1.feb'))
A = np.array([f.nodes[23700 + k] for k in range(1, 14)])       # arcus nodes 1-13 (left)
B = np.array([f.nodes[n] for n in PVW])
L = np.linalg.norm(B - A, axis=1)
def spacing(X):
    d = np.linalg.norm(np.diff(X, axis=0), axis=1)
    return np.r_[d[0], (d[:-1] + d[1:]) / 2, d[-1]] if len(d) > 1 else d
w = (spacing(A) + spacing(B)) / 2
print('connector lengths L [mm]:', np.round(L, 1))
print('strip widths w [mm]:', np.round(w, 1), ' total', round(w.sum(), 1))
u = np.linspace(2, 40, 39)
target = np.array([np.interp(x, uk, Fk) * 13 for x in u])       # all 13 connectors (same table), per side
def fan_force(c, m, x):
    lam = 1 + x / L
    return np.sum((c / m) * (lam ** (m - 1) - lam ** (-m / 2 - 1)) * w * t)
best = None
for m in np.arange(1.0, 8.01, 0.05):
    g = np.array([fan_force(1.0, m, x) for x in u])                # force for c = 1 (linear in c)
    c = np.sum(target * g / target ** 2) / np.sum(g * g / target ** 2)
    err = np.sqrt(np.mean(((c * g - target) / target) ** 2))
    if best is None or err < best[0]:
        best = (err, c, m)
err, c, m = best
mu0 = c / 2
print(f'best 1-term Ogden: c1 = {c:.5f} MPa, m1 = {m:.2f} (mu0 = {mu0:.5f} MPa, E0 ~ {3 * mu0:.4f} MPa), rms rel. error {100 * err:.1f} %')
for x in (4.7, 10, 20, 30, 40):
    lam = 1 + x / L
    fan = np.sum((c / m) * (lam ** (m - 1) - lam ** (-m / 2 - 1)) * w * t)
    print(f'  u = {x:4.1f} mm: connectors {13 * np.interp(x, uk, Fk):6.2f} N, soft fan {fan:6.2f} N, '
          f'current fan (E 21, linear strip) ~{np.sum(21 * (x / L) * w * t):7.1f} N')
