"""Batch 337 (2026-10-19; the user: option 3, "test the fibre material on a small solid strip, then build the solid fibre strap
on L288c (no contact yet)"). The curve-tip -> PM fan as a SOLID FIBRE STRAP: from SRC (= L288c + the fan as chains in a sheet,
build_batch334.py), (1) the fan chains' directions at rest (PM end -> curve end) are read and the chains removed; (2) the two
fan sheets become one penta6 layer TH thick (build_batch320's routine, no contact; the strap keeps the sheet's nodes, so its
curve edge stays on the spring-patch anchors); (3) the strap material: uncoupled solid mixture, Mooney-Rivlin matrix C1 (c2 0)
+ a tension-only fibre-pow-linear-uncoupled (E FE, beta 2, lam0 1.01) along the chains' direction per element (MeshData
mat_axis, a = the direction of the two nearest chains, projected onto the element's plane), mixture k K.
Solid strip test (shell_tests/fibre_solid, penta6 10 x 4 x 0.49 mm, 2026-10-19): with matrix c1 0.00833 and k 10, fibre E 27.5
pulls 3.88 at 10 % strain (E 22 -> 3.29; ~0.76 of the plain E A strain estimate: the volume change at k 10), compression
~0.006 (650x less), across the fibres 0.011; k >= 100 fails in tension; the coupled mixture failed at start-up. The springs'
strip fit (u 0.25-4 mm, the 5 chains a side) gives a linear fibre E 22.2-22.4 MPa (rms 5 %) -> input E 22.3 / 0.76 = 29.5.
NOT IN SOURCE (a strap in place of connectors that are themselves not in the source).
usage: py -3.10 build_batch337.py SRC_RUN NEW_RUN [--th 0.49] [--fe 29.5] [--c1 0.00833] [--k 10]"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_batch111 import Model111  # noqa: E402
from build_batch332 import chains  # noqa: E402
from variants2 import renumber_discrete_materials  # noqa: E402
from variants import emit  # noqa: E402
import build_batch320  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
SIDES = ('left', 'right')
# family -> [(sheet, [chain sets])] (2026-10-19 L342: generalised from the fan to any chains-in-a-sheet family)
FAMS = {'fan': [(f'PM_avw_bottom_{s}_sheet', [f'curve_end_PM_conn_PM_avw_bottom_{s}_chain', f'PM_avw_bottom_{s}_conn_PM_avw_bottom_{s}_chain'])
                for s in SIDES],
        'AVW-Para': [(f'AVW-Para-{s}_sheet', [f'AVW-Para-{s}_conn_chain']) for s in 'LR']}
SOLID_FACTOR = 0.76   # strip test: the fibre pulls ~0.76 of E A strain at k 10 (volume change)


def fit_fibre(con, area, u0, u1, th):
    """Linear fibre E (MPa) whose strap pulls like the chains at uniform end displacement u (strip model as loft_fit_all:
    each chain a strip of its rest length L_i, width area / sum L, thickness th); returns (E_linear, rms rel error)."""
    u = np.linspace(u0, u1, 40)
    w = area / sum(c[0] for c in con)
    tg = sum(sc * np.interp(u / n, pts[:, 0], pts[:, 1]) for L, n, pts, sc in con)
    g = sum((u / L - 0.005) * w * th for L, *_ in con)
    E = np.sum(g / tg) / np.sum((g / tg) ** 2)
    return E, float(np.sqrt(np.mean((E * g / tg - 1) ** 2)))


def main(src, new, th=0.49, fe=29.5, c1=0.00833, k=10.0, fam='fan', mid=None, u0=0.25, u1=4.0):
    specs = FAMS[fam]
    # (1) chain directions, then remove the chains (intermediate NEW + 'a', built only)
    m = Model111(os.path.join(RUNS, src, src + '.feb'))
    X = m.nodes()
    rails = {sh: [] for sh, _ in specs}
    fes = {}
    disc = m.root.find('Discrete')
    for sh, sets in specs:
        con = []
        for s in sets:
            pairs, mat, b, ds = m._spring_set(s)
            pts = np.array([[float(v) for v in p_.text.split(',')] for p_ in mat.find('force').find('points')])
            sc = float(mat.findtext('scale'))
            for c in chains(pairs):
                A, B = X[c[0]], X[c[-1]]
                rails[sh].append((np.array([X[n] for n in c]), (B - A) / np.linalg.norm(B - A)))
                con.append((sum(np.linalg.norm(X[c[q]] - X[c[q + 1]]) for q in range(len(c) - 1)), len(c) - 1, pts, sc))
            disc.remove(b); disc.remove(mat); m.mesh.remove(ds)
            m.log.append(f'{s}: its {len(pairs)} chain springs removed')
        if fe == 'fit':
            tris = next(e for e in m.mesh.findall('Elements') if e.get('name') == sh)
            area = sum(0.5 * np.linalg.norm(np.cross(X[q[1]] - X[q[0]], X[q[2]] - X[q[0]]))
                       for q in ([int(v) for v in e.text.split(',')] for e in tris))
            El, err = fit_fibre(con, area, u0, u1, th)
            fes[sh] = El / SOLID_FACTOR
            m.log.append(f'{sh}: fibre fit to its {len(con)} chains over u {u0:g}-{u1:g} mm: linear E {El:.4g} MPa (rms {err * 100:.0f} %), '
                         f'input E {fes[sh]:.4g} (/ {SOLID_FACTOR})')
        else:
            fes[sh] = float(fe)
    renumber_discrete_materials(m)
    mid = mid or new.split('_')[0] + 'a_' + '_'.join(new.split('_')[1:])
    assert not os.path.exists(os.path.join(RUNS, mid)), mid
    emit(mid, m)
    # (2) the fan sheets as one penta6 layer (build_batch320, no contact)
    sheets = tuple(sh for sh, _ in specs)
    build_batch320.main(mid, new, th, False, 1.0, 3.0, sheets, ('none',), 'fanstrap' if fam == 'fan' else fam.replace('-', '') + 'strap')
    # (3) the strap material + fibre directions
    path = os.path.join(RUNS, new, new + '.feb')
    t = open(path, encoding='utf-8').read()
    nodes = {int(i): np.array([float(v) for v in c.split(',')]) for i, c in re.findall(r'<node id="(\d+)">([^<]+)</node>', t)}
    data, log = [], []
    for sh, _ in specs:
        dom = f'{sh}_solid'
        body = re.search(rf'<Elements type="penta6" name="{dom}">(.*?)</Elements>', t, re.S).group(1)
        els = [[int(v) for v in c.split(',')] for c in re.findall(r'<elem id="\d+">([^<]+)</elem>', body)]
        rows, ang = [], []
        for e in els:
            P = np.array([nodes[v] for v in e[:3]]); cen = P.mean(0)
            nrm = np.cross(P[1] - P[0], P[2] - P[0]); nrm /= np.linalg.norm(nrm)
            dist = [np.min(np.linalg.norm(R - cen, axis=1)) for R, _ in rails[sh]]
            j = np.argsort(dist)[:2]
            w = 1.0 / (np.array(dist)[j] + 1e-6)
            a = sum(wi * rails[sh][ji][1] for wi, ji in zip(w, j)); a = a - np.dot(a, nrm) * nrm; a /= np.linalg.norm(a)
            d = np.cross(nrm, a)
            ang.append(np.degrees(np.arccos(abs(np.dot(a, rails[sh][j[0]][1])))))
            rows.append(f'\t\t\t<elem lid="{len(rows) + 1}"><a>{a[0]:.6f},{a[1]:.6f},{a[2]:.6f}</a><d>{d[0]:.6f},{d[1]:.6f},{d[2]:.6f}</d></elem>\n')
        data.append(f'\t\t<ElementData type="mat_axis" elem_set="{dom}">\n' + ''.join(rows) + '\t\t</ElementData>\n')
        log.append(f'{dom}: {len(els)} penta6, fibre direction per element = the 2 nearest chains (projected onto the element, '
                   f'{max(ang):.1f} deg at most from the nearest chain)')
        mname = re.search(rf'<SolidDomain name="{dom}" mat="([^"]+)"', t).group(1)
        old = re.search(rf'<material id="(\d+)" name="{mname}" type="[^"]+">.*?</material>', t, re.S)
        newmat = (f'<material id="{old.group(1)}" name="{mname}" type="uncoupled solid mixture">\n\t\t\t<density>1.06e-09</density>\n'
                  f'\t\t\t<k>{k:g}</k>\n\t\t\t<solid type="Mooney-Rivlin">\n\t\t\t\t<c1>{c1:g}</c1>\n\t\t\t\t<c2>0</c2>\n\t\t\t</solid>\n'
                  f'\t\t\t<solid type="fiber-pow-linear-uncoupled">\n\t\t\t\t<fiber type="vector">1,0,0</fiber>\n\t\t\t\t<E>{fes[sh]:.4g}</E>\n'
                  f'\t\t\t\t<beta>2</beta>\n\t\t\t\t<lam0>1.01</lam0>\n\t\t\t</solid>\n\t\t</material>')
        t = t.replace(old.group(0), newmat, 1)
    assert '<MeshData>' not in t or True
    if '<MeshData>' in t:
        i = t.index('</MeshData>'); t = t[:i] + ''.join(data) + '\t' + t[i:]
    else:
        i = t.index('</MeshDomains>') + len('</MeshDomains>\n'); t = t[:i] + '\t<MeshData>\n' + ''.join(data) + '\t</MeshData>\n' + t[i:]
    open(path, 'w', encoding='utf-8').write(t)
    msg = (f'NOT IN SOURCE: the {fam} family as a SOLID FIBRE STRAP (build_batch337.py): chains removed ({src} -> {mid}, built only), '
           f'sheets as one penta6 layer {th:g} mm (build_batch320), material uncoupled solid mixture k {k:g}, Mooney-Rivlin c1 {c1:g}, '
           f'tension-only fiber-pow-linear-uncoupled E {", ".join(f"{k_} {v:.4g}" for k_, v in fes.items())} beta 2 lam0 1.01; ' + '; '.join(log))
    with open(path + '.changes.txt', 'a') as f:
        f.write(msg + '\n')
    print('\n'.join(m.log)); print(msg)


if __name__ == '__main__':
    a = sys.argv[1:]
    g = lambda key, d: float(a[a.index(key) + 1]) if key in a else d   # noqa: E731
    fe = a[a.index('--fe') + 1] if '--fe' in a else '29.5'
    main(a[0], a[1], g('--th', 0.49), fe if fe == 'fit' else float(fe), g('--c1', 0.00833), g('--k', 10.0),
         a[a.index('--fam') + 1] if '--fam' in a else 'fan', a[a.index('--mid') + 1] if '--mid' in a else None, g('--u0', 0.25), g('--u1', 4.0))
