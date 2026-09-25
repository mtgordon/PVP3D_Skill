"""Task 3 checks. (a) Abaqus BC-VW-mid (XSYMM on _PickedSet1359) vs FEBio BC-VW-mid / StabilizationAnchors.
(b) Abaqus CL/USL *Rigid Body (beam fully fixed via ref node) vs FEBio fan nodes on the beam line."""
import sys
import numpy as np
from abq_surf import lines
from febmodel import Feb
from inpmap import part_nodes

feb = Feb(sys.argv[1])
X = feb.nodes

def nset_inst(name):
    out = []; grab = gen = False
    for ln in lines():
        s = ln.strip()
        if s.startswith('*'):
            if grab: break
            if s.lower().startswith('*nset') and ('nset=%s,' % name in s or s.endswith('nset=%s' % name)):
                grab = True; gen = 'generate' in s.lower()
            continue
        if grab and s:
            v = [int(x) for x in s.replace(' ', '').split(',') if x]
            out += list(range(v[0], v[1] + 1, v[2] if len(v) > 2 else 1)) if gen else v
    return out

vw = part_nodes('VW-PeB')
ids = nset_inst('_PickedSet1359')
P = np.array([vw[i] for i in ids])
print(f'(a) Abaqus _PickedSet1359: {len(ids)} VW-PeB nodes; x range {P[:,0].min():.4f}..{P[:,0].max():.4f}')
fid = np.array(sorted(X)); FX = np.array([X[i] for i in fid])
def match(p):
    d = np.linalg.norm(FX - p, axis=1); k = int(np.argmin(d)); return int(fid[k]), float(d[k])
m = [match(p) for p in P]
print('    max match distance %.2e' % max(d for _, d in m))
abq_mid = {i for i, _ in m}
vwmid = set(feb.nodesets['BC-VW-mid']); anch = set(feb.nodesets['StabilizationAnchorsSet'])
print(f'    FEBio BC-VW-mid {len(vwmid)} nodes; == Abaqus set: {vwmid == abq_mid}; overlap {len(vwmid & abq_mid)}')
print(f'    StabilizationAnchors {len(anch)}; overlap with VW-mid {len(anch & abq_mid)}')
# which nodes do the stab springs connect (spring node -> anchor)?
sp = [(a, b) for n, pairs in feb.discsets.items() if n.startswith('stab_') for a, b in pairs]
live = {a if b in anch else b for a, b in sp}
print(f'    stab springs {len(sp)}; spring nodes on the Abaqus VW-mid set: {len(live & abq_mid)}; '
      f'spring nodes elsewhere: {len(live - abq_mid)}')
A = np.array([X[a] for a in live]); print('    spring-node x range %.3f..%.3f' % (A[:,0].min(), A[:,0].max()))
used = feb.nodes_in_elements()
print(f'    Abaqus VW-mid nodes used by FEBio elements: {len(abq_mid & used)}')
dom = {}
for name, (et, d) in feb.elem_blocks.items():
    for c in d.values():
        for n in c:
            if n in abq_mid: dom.setdefault(name, set()).add(n)
print('    by domain:', {k: len(v) for k, v in dom.items()})
# (b) CL/USL beams
fixed = set()
for n, s in feb.nodesets.items():
    if n.startswith(('BC-CL', 'BC-USL')): fixed |= set(s)
for part, fan in (('CL_Left', 'CL-L_fan'), ('CL_Right', 'CL-R_fan'), ('USL_Left', 'USL-L_fan'), ('USL_Right', 'USL-R_fan')):
    B = part_nodes(part); bn = np.array([B[i] for i in sorted(B)])
    fn = feb.domain_nodes(fan)
    def dline(p):
        best = 1e9
        for a, b in zip(bn[:-1], bn[1:]):
            u = b - a; t = np.clip((p - a) @ u / (u @ u), 0, 1); best = min(best, np.linalg.norm(p - a - t * u))
        return best
    near = [(n, dline(X[n])) for n in fn]
    on = [n for n, d in near if d < 0.6]
    print(f'(b) {fan}: {len(fn)} nodes, {len(on)} within 0.6 mm of the {part} beam line, fixed: '
          f'{sum(n in fixed for n in on)}; nearest unfixed: {min((d for n, d in near if n not in fixed), default=None):.3f} mm')
