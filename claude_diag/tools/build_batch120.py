"""Batch 120 (2026-09-29, the user's request): the levator (LA) and perineal body (PeB) with Jing, Ashton-Miller &
DeLancey 2012 (J Biomech 45:455-460; Paper with Levator Tissue Data (Dejun Jing).pdf;
text claude_diag/scratch_2026-09-29/jing_levator.txt, eq. 1 read from the rendered page 3):
  W = C (I1 - 3) + k1 / (2 k2) [exp(k2 (lam^2 - 1)^2) - 1]  (fibre term for lam >= 1 only; the paper prints (I1 - 1)
  and no "- 1" in the bracket, taken as the standard Holzapfel form), quasi-linear viscoelastic with long-term G_inf;
  Table 1 (fitted to biaxial tests, scaled for pregnancy): LA C 0.181, k1 0.083 MPa, k2 0.32, G_inf 0.22;
  PeB C 0.293, k1 0.037, k2 0.65, G_inf 0.17.
Our loads are slow and held, so the RELAXED law is the comparison: C and k1 x G_inf (QLV scales the whole elastic
response): LA C 0.03982, k1 0.01826 MPa; PeB C 0.04981, k1 0.00629 MPa (k2 unchanged).
In FEBio: an uncoupled solid mixture of Mooney-Rivlin (c1 = C, c2 = 0: the neo-Hookean matrix) and
fiber-exp-pow-uncoupled (ksi = k1, alpha = k2, beta = 2: the same fibre energy), k = 43 x the matrix shear modulus 2C
(the lines' LA ratio; PeB: its current k 5.7, the same small-strain modulus), density 1.06e-9.
Fibre directions (the paper follows Shobeiri et al. 2008; our LA mesh has none, NOT IN SOURCE):
  LA_PCMPRM, LA_PCM (the pubovisceral part): along the hiatus's free edge (the tangent of the nearest free-edge segment
  of the PCM / PRM, projected into the element);
  LA_ICM, LA_ICM_tri (iliococcygeus): toward the nearest midline-raphe node (BC-LA-mid), projected into the element;
  PeB: transverse (global x).
Runs, one change each from the springs line's base L91_newline_vwyeoh_rhoi0_seamspr001:
  L120_jing_la          the LA with the relaxed Jing law (fibres as above)
  L120_jing_peb         the PeB with the relaxed Jing law
  L120_jing_la_peb      both
  L120_jing_la_inst     the LA with the instantaneous (unrelaxed) Jing law: the upper bracket
usage: py -3.10 build_batch120.py [NAME ...]"""
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import Model7, JOBS, emit  # noqa: E402
from febmodel import Feb  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
BASE = 'L91_newline_vwyeoh_rhoi0_seamspr001'
WANT = set(sys.argv[1:])
LA_DOMS = ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri')
JING = {'LA': dict(C=0.181, k1=0.083, k2=0.32, G=0.22), 'PeB': dict(C=0.293, k1=0.037, k2=0.65, G=0.17)}
LABEL = 'NOT IN SOURCE (Jing et al. 2012 tissue law)'


def mixture(name, mid, C, k1, k2, k, fibre_vec=None):
    m = ET.Element('material', id=str(mid), name=name, type='uncoupled solid mixture')
    ET.SubElement(m, 'density').text = '1.06e-09'
    ET.SubElement(m, 'k').text = f'{k:.6g}'
    s = ET.SubElement(m, 'solid', type='Mooney-Rivlin')
    ET.SubElement(s, 'c1').text = f'{C:.6g}'
    ET.SubElement(s, 'c2').text = '0'
    s = ET.SubElement(m, 'solid', type='fiber-exp-pow-uncoupled')
    ET.SubElement(s, 'fiber', type='vector').text = fibre_vec or '1,0,0'
    ET.SubElement(s, 'ksi').text = f'{k1:.6g}'
    ET.SubElement(s, 'alpha').text = f'{k2:.6g}'
    ET.SubElement(s, 'beta').text = '2'
    return m


def seg_tangent(p, segs, X):
    best, t = 1e30, None
    for a, b in segs:
        A, B = X[a], X[b]
        ab = B - A
        s = np.clip(np.dot(p - A, ab) / np.dot(ab, ab), 0, 1)
        d = np.linalg.norm(p - (A + s * ab))
        if d < best:
            best, t = d, ab / np.linalg.norm(ab)
    return t


def la_fibres(f):
    """{domain: [(a, d) per element in file order]}: the fibre (a) and a second axis (d, the element normal)."""
    X = f.nodes
    # the LA's free edges: boundary edges of all LA elements whose nodes are not fixed, tied or on another domain
    ec = Counter()
    for dom in LA_DOMS:
        for c in f.elem_blocks[dom][1].values():
            for i in range(len(c)):
                ec[tuple(sorted((c[i], c[(i + 1) % len(c)])))] += 1
    held = {n for b in f.root.find('Boundary') for n in f.nodesets.get(b.get('node_set'), [])}
    held |= {int(n.get('id')) for c in f.root.find('Constraints') for lc in c.findall('linear_constraint') for n in lc.findall('node')}
    pv = {n for dom in ('LA_PCMPRM', 'LA_PCM') for c in f.elem_blocks[dom][1].values() for n in c}
    free_pv = [e for e, k in ec.items() if k == 1 and e[0] in pv and e[1] in pv and e[0] not in held and e[1] not in held]
    raphe = np.array([X[n] for n in f.nodesets['BC-LA-mid']])
    out = {}
    for dom in LA_DOMS:
        rows = []
        for c in f.elem_blocks[dom][1].values():
            P = np.array([X[n] for n in c])
            cen = P.mean(axis=0)
            nrm = np.cross(P[1] - P[0], P[2] - P[0]) if len(c) == 3 else np.cross(P[2] - P[0], P[3] - P[1])
            nrm /= np.linalg.norm(nrm)
            if dom in ('LA_PCMPRM', 'LA_PCM'):
                v = seg_tangent(cen, free_pv, X)
            else:
                v = raphe[int(np.argmin(np.linalg.norm(raphe - cen, axis=1)))] - cen
            v = v - (v @ nrm) * nrm
            v /= np.linalg.norm(v)
            rows.append((v, nrm))
        out[dom] = rows
    return out, len(free_pv)


class Model120(Model7):
    def jing_la(self, relaxed=True):
        p = JING['LA']
        g = p['G'] if relaxed else 1.0
        C, k1 = p['C'] * g, p['k1'] * g
        mats = self.root.find('Material')
        mid = max(int(m.get('id')) for m in mats) + 1
        name = 'LA_Jing2012' + ('_relaxed' if relaxed else '_instant')
        mats.append(mixture(name, mid, C, k1, p['k2'], 43 * 2 * C))
        for d in self.root.find('MeshDomains'):
            if d.get('name') in LA_DOMS:
                d.set('mat', name)
        f = Feb(self.src)
        fib, nfree = la_fibres(f)
        md = self.root.find('MeshData')
        if md is None:
            md = ET.Element('MeshData')
            kids = list(self.root)
            self.root.insert(kids.index(self.root.find('MeshDomains')) + 1, md)
        for dom, rows in fib.items():
            ed = ET.SubElement(md, 'ElementData', type='mat_axis', elem_set=dom)
            for i, (a, d) in enumerate(rows, 1):
                e = ET.SubElement(ed, 'elem', lid=str(i))
                ET.SubElement(e, 'a').text = ','.join(f'{v:.6f}' for v in a)
                ET.SubElement(e, 'd').text = ','.join(f'{v:.6f}' for v in d)
        self.log.append(f'{LABEL}: the LA (LA_PCMPRM, LA_PCM, LA_ICM, LA_ICM_tri) from LA_Yamada50pct_Yeoh (the healthy '
                        f'Yamada 100 % Yeoh) to {name}: uncoupled solid mixture, Mooney-Rivlin c1 {C:.5g} c2 0 + '
                        f'fiber-exp-pow-uncoupled ksi {k1:.5g} alpha {p["k2"]} beta 2, k {43 * 2 * C:.4g}, density 1.06e-9 '
                        f'(Jing 2012 Table 1 {"x G_inf " + str(p["G"]) + ", the relaxed law" if relaxed else "as given, instantaneous"}); '
                        f'fibres per element (mat_axis): PCM / PRM along the hiatus free edge ({nfree} free edges), ICM '
                        f'toward the nearest BC-LA-mid (raphe) node, in the element plane')

    def jing_peb(self, relaxed=True):
        p = JING['PeB']
        g = p['G'] if relaxed else 1.0
        C, k1 = p['C'] * g, p['k1'] * g
        mats = self.root.find('Material')
        old = next(m for m in mats if m.get('name') == 'PeB-Vagina500%stiffer')
        k_old = float(old.findtext('k'))
        mid = max(int(m.get('id')) for m in mats) + 1
        name = 'PeB_Jing2012' + ('_relaxed' if relaxed else '_instant')
        mats.append(mixture(name, mid, C, k1, p['k2'], k_old, '1,0,0'))
        for d in self.root.find('MeshDomains'):
            if d.get('mat') == 'PeB-Vagina500%stiffer':
                d.set('mat', name)
        self.log.append(f'{LABEL}: the perineal body (_PickedSet66) from PeB-Vagina500%stiffer (the Yeoh refit) to {name}: '
                        f'uncoupled solid mixture, Mooney-Rivlin c1 {C:.5g} c2 0 + fiber-exp-pow-uncoupled ksi {k1:.5g} '
                        f'alpha {p["k2"]} beta 2 along global x (transverse), k {k_old:g} (unchanged), density 1.06e-9')


BUILDS = [('L120_jing_la', lambda m: m.jing_la()),
          ('L120_jing_peb', lambda m: m.jing_peb()),
          ('L120_jing_la_peb', lambda m: (m.jing_la(), m.jing_peb())),
          ('L120_jing_la_inst', lambda m: m.jing_la(relaxed=False))]

if __name__ == '__main__':
    for name, fn in BUILDS:
        if WANT and name not in WANT:
            continue
        if os.path.exists(os.path.join(RUNS, name)):
            print(f'{name} exists; not overwriting')
            continue
        m = Model120(os.path.join(RUNS, BASE, BASE + '.feb'))
        fn(m)
        emit(name, m)
