"""The deepest curve-through-strap crossings with node / face / tie details (2026-10-23).
usage: py -3.10 strap_first_poke.py RUN T1,T2"""
import sys, os, re, numpy as np
from collections import Counter, defaultdict
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
dn=defaultdict(set)
for k,(n,c) in E.items():
    for v in c: dn[v].add(n)
sec={frozenset(map(int,f.split(','))) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="curve_strap_secondary">(.*?)</Surface>',t,re.S).group(1))}
prim={int(v) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="curve_strap_primary">(.*?)</Surface>',t,re.S).group(1)) for v in f.split(',')}
cons=set()
for m in re.finditer(r'<constraint name="([^"]+)"(.*?)</constraint>',t,re.S):
    cons|={int(a) for a in re.findall(r'<linear_constraint>\s*<node id="(\d+)" bc="x"',m.group(2))}
curve={v for v,s in dn.items() if any(q.startswith('canal_wrap') for q in s)}
strapn={v for v,s in dn.items() if any(q.startswith('PM_avw_bottom') for q in s)}
st,se=strap_surface(x); ct,ce=surface(x)
X0=x.X; b={c[:3] for c in crossings(X0,st,ce)}
tt=np.array([x.states[i][0] for i in idx])
for T in [float(v) for v in sys.argv[2].split(',')]:
    i=idx[int(np.argmin(abs(tt-T)))]; P=X0+x.var(i,'displacement')
    cr=[c for c in crossings(P,st,ce) if c[:3] not in b]
    print('t %.4f'%x.states[i][0], len(cr),'curve edges through strap faces')
    for a,bb,j,dep in sorted(cr,key=lambda c:-c[3])[:6]:
        A,B=int(mp[a]),int(mp[bb]); tri=[int(mp[q]) for q in st[j]]
        face=next((f for f in sec if set(tri)<=f),None)
        shared=[q for q in tri if q in curve]
        print(f'  edge {A}-{B} depth {dep:.2f}  rest {np.round((X0[a]+X0[bb])/2,1)}  in prim {A in prim}/{B in prim}  tied {A in cons}/{B in cons} | strap tri {tri} in curve_strap_secondary {face is not None}; tri nodes that are CURVE nodes (tied rim) {shared}; strap-tied {[q in cons for q in tri]}')
