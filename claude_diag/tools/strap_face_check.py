"""Which strap faces (top / bottom / side) the curve pierces and whether they and the nodes are in curve_strap (2026-10-23).
usage: py -3.10 strap_face_check.py RUN T1,T2"""
import sys, os, re, numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from strap_cross import strap_surface, surface, crossings, RUNS, PF
from xplt_reader import Xplt, converged_state_indices
from junction_view_data import load
from scipy.spatial import cKDTree
run=sys.argv[1]
x=Xplt(os.path.join(RUNS,run,run+'.xplt')); idx=converged_state_indices(x,os.path.join(RUNS,run,run+'.log'))
t,FX,E=load(run); ids=list(FX); d,ii=cKDTree(np.array([FX[k] for k in ids])).query(x.X); mp=np.array([ids[k] for k in ii])
sec={frozenset(map(int,f.split(','))) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="curve_strap_secondary">(.*?)</Surface>',t,re.S).group(1))}
fsp={frozenset(map(int,f.split(','))) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="fanstrap_curves_primary">(.*?)</Surface>',t,re.S).group(1))}
prim={int(v) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="curve_strap_primary">(.*?)</Surface>',t,re.S).group(1)) for v in f.split(',')}
W=np.vstack([dd['conn'][:,:6] for dd in x.domains if (dd.get('name') or '').startswith('PM_avw_bottom_')])
ftype={}
for w in W:
    for k,f in enumerate(PF):
        ftype[frozenset(mp[w[list(f)]])]=('bottom','top','side','side','side')[k]
st,se=strap_surface(x); ct,ce=surface(x)
def which(tri):
    s=frozenset(mp[tri])
    for f,ty in ftype.items():
        if s<=f: return ty, any(s<=g for g in sec)
tt=np.array([x.states[i][0] for i in idx]); X0=x.X
b={c[:3] for c in crossings(X0,st,ce)}
print('curve_strap_secondary', len(sec), 'facets; types', Counter(ftype.get(f,'?') for f in sec), '; curve_strap_primary nodes', len(prim))
for T in sys.argv[2].split(','):
    i=idx[int(np.argmin(abs(tt-float(T))))]; P=X0+x.var(i,'displacement')
    cr=[c for c in crossings(P,st,ce) if c[:3] not in b]
    print('t %.3f'%x.states[i][0], 'pierced strap faces (type, in curve_strap_secondary):', Counter(which(st[j]) for a,bb,j,dd in cr),
          '; piercing curve edge nodes in curve_strap_primary:', Counter((int(mp[a]) in prim, int(mp[bb]) in prim) for a,bb,j,dd in cr))
