"""Which Yeoh term controls each vaginal wall (2026-10-07, the user: "do the check to see which parameter is controlling each
wall"). FEBio's Yeoh is uncoupled: W = c1 (I1~ - 3) + c2 (I1~ - 3)^2 + U(J), I1~ = J^(-2/3) tr(F F^T). Its stiffness against
I1~ is dW/dI1~ = c1 + 2 c2 (I1~ - 3): the c1 term is constant, the c2 term grows with stretch. Per wall element at its
centroid (hex8 shape-function gradient of the plot's nodal displacements) at the chosen state: I1~ - 3, the equivalent
uniaxial strain (lambda - 1 with lambda^2 + 2 / lambda = I1~), and the c2 share 2 c2 (I1~-3) / (c1 + 2 c2 (I1~-3)).
Per wall (AVW, cervix, PVW, + the lips if present) and per canal third (distal = nearest the hymen): median / p90 / max.
usage: py -3.10 wall_stretch.py RUN [--t 1]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.expanduser('~/.claude/skills/abaqus-febio-fea-pipeline/scripts'))
from xplt_reader import Xplt  # noqa: E402
from variants7 import JOBS  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
WALLS = {'_PickedSet347': 'AVW', '_PickedSet346': 'cervix', '_PickedSet64': 'PVW',
         'AVW_ext': 'AVW lip', 'cervix_ext': 'cervix lip', 'PVW_ext': 'PVW lip'}
DN = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]) / 8.0


def lam_of(i1):
    """uniaxial stretch with lambda^2 + 2/lambda = I1 (lambda >= 1 branch)"""
    lo, hi = 1.0, 5.0
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if m * m + 2 / m < i1 else (lo, m)
    return (lo + hi) / 2


if __name__ == '__main__':
    args = sys.argv[1:]
    tw = 1.0
    if '--t' in args:
        i = args.index('--t'); tw = float(args[i + 1]); args = args[:i] + args[i + 2:]
    run = args[0]
    t = open(os.path.join(RUNS, run, run + '.feb'), encoding='utf-8').read()
    X0 = {int(m.group(1)): np.array([float(v) for v in m.group(2).split(',')])
          for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t)}
    els = {}
    for nm in WALLS:
        m = re.search(r'<Elements type="hex8" name="%s">(.*?)</Elements>' % re.escape(nm), t, re.S)
        if m:
            els[nm] = np.array([[int(v) for v in e.split(',')] for e in re.findall(r'>([\d,]+)</elem>', m.group(1))])
    mat = {}
    for nm in els:
        dom = re.search(r'<SolidDomain name="%s" mat="([^"]+)"' % re.escape(nm), t).group(1)
        mb = re.search(r'<material id="\d+" name="%s" type="Yeoh">(.*?)</material>' % re.escape(dom), t, re.S).group(1)
        mat[nm] = tuple(float(re.search(r'<%s>([^<]+)</%s>' % (k, k), mb).group(1)) for k in ('c1', 'c2'))
    x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
    si = min(range(len(x.states)), key=lambda i: abs(x.states[i][0] - tw))
    u = np.asarray(x.var(si, 'displacement'), dtype=float)
    idx = {int(n): i for i, n in enumerate(x.node_ids)}
    # canal axis and the hymen end: the AVW's principal axis, distal = the end nearer the PM plane nodes is not known here,
    # so use the end with the larger mean x-height drop... simpler: distal = the end of the axis nearest the PeB (_PickedSet66)
    avw = np.array([X0[v] for v in els['_PickedSet347'].ravel()])
    c0 = avw.mean(0)
    ax = np.linalg.svd(avw - c0, full_matrices=False)[2][0]
    peb = re.search(r'<Elements type="hex8" name="_PickedSet66">(.*?)</Elements>', t, re.S).group(1)
    pc = np.mean([X0[int(v)] for e in re.findall(r'>([\d,]+)</elem>', peb) for v in e.split(',')], axis=0)
    if (pc - c0) @ ax < 0:
        ax = -ax                       # + along ax = towards the hymen
    allp = np.array([X0[v] for e in els.values() for v in e.ravel()])
    s_all = (allp - c0) @ ax
    lo, hi = np.percentile(s_all, [1, 99])
    print(f'{run} at t {x.states[si][0]:.3f}: Yeoh dW/dI1~ = c1 + 2 c2 (I1~-3); "c2 share" = 2 c2 (I1~-3) / that')
    print('| wall | region | n | I1~-3 median / p90 / max | uniaxial-equivalent strain median / p90 / max | c2 share median / p90 |')
    print('|---|---|---|---|---|---|')
    for nm, E in els.items():
        c1, c2 = mat[nm]
        P0 = np.array([[X0[v] for v in e] for e in E])
        Uu = np.array([[u[idx[v]] for v in e] for e in E])
        J0 = np.einsum('eai,aj->eij', P0, DN)              # dX/dxi
        G = np.einsum('eai,aj->eij', Uu, DN)               # du/dxi
        F = np.eye(3) + np.einsum('eij,ejk->eik', G, np.linalg.inv(J0))
        J = np.linalg.det(F)
        I1 = np.einsum('eij,eij->e', F, F) * np.abs(J) ** (-2 / 3)
        d = np.maximum(I1 - 3, 0)
        strain = np.array([lam_of(3 + v) - 1 for v in d])
        share = 2 * c2 * d / (c1 + 2 * c2 * d)
        s = (P0.mean(1) - c0) @ ax
        thirds = np.digitize(s, [lo + (hi - lo) / 3, lo + 2 * (hi - lo) / 3])
        for reg, mask in (('all', np.ones(len(E), bool)), ('proximal', thirds == 0), ('middle', thirds == 1),
                          ('distal', thirds == 2)):
            if mask.sum() == 0:
                continue
            q = lambda a: np.percentile(a[mask], [50, 90])  # noqa: E731
            print(f'| {WALLS[nm]} | {reg} | {mask.sum()} | {q(d)[0]:.3f} / {q(d)[1]:.3f} / {d[mask].max():.3f} | '
                  f'{q(strain)[0] * 100:.0f} % / {q(strain)[1] * 100:.0f} % / {strain[mask].max() * 100:.0f} % | '
                  f'{q(share)[0] * 100:.0f} % / {q(share)[1] * 100:.0f} % |')
    print(f'(c1 term = c2 term at I1~-3 = c1 / (2 c2) = {mat["_PickedSet347"][0] / (2 * mat["_PickedSet347"][1]):.3f}, '
          f'a uniaxial-equivalent strain of {(lam_of(3 + mat["_PickedSet347"][0] / (2 * mat["_PickedSet347"][1])) - 1) * 100:.0f} %)')
