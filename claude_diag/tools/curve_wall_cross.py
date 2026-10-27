"""Curve skin edges (exterior faces of the canal_wrap_hex domains) crossing WALL-ONLY faces (walls_LA_primary faces that are not all curve-hex nodes) at the last state.
usage: py -3.10 curve_wall_cross.py RUN (2026-10-26, the colleague)"""
import os,sys,re,numpy as np
sys.path.insert(0,'/home/claude/pvp3d_skill/claude_diag/tools')
from curve_cross import surface, crossings, RUNS
from paths import SKILL_SCRIPTS; sys.path.insert(0,SKILL_SCRIPTS)
from xplt_reader import Xplt, converged_state_indices
run=sys.argv[1]
x=Xplt(os.path.join(RUNS,run,run+'.xplt')); idx=converged_state_indices(x,os.path.join(RUNS,run,run+'.log'))
tris,edges=surface(x); X0=x.X
t=open(os.path.join(RUNS,run,run+'.feb')).read()
ids=[int(v) for v in re.findall(r'<node id="(\d+)"',re.search(r'<Mesh>(.*?)</Mesh>',t,re.S).group(1))]; ix={n:i for i,n in enumerate(ids)}
H=np.vstack([d['conn'][:,:8] for d in x.domains if (d.get('name') or '').startswith('canal_wrap_hex')])
cnodes=set(int(v) for v in H.ravel())
W=[[ix[int(v)] for v in f.split(',')] for f in re.findall(r'>([^<]+)</(?:quad4|tri3)>',re.search(r'<Surface name="walls_LA_primary">(.*?)</Surface>',t,re.S).group(1))]
wonly=[f for f in W if not all(n in cnodes for n in f)]
wt=np.array([[f[0],f[1],f[2]] for f in wonly]+[[f[0],f[2],f[3]] for f in wonly if len(f)==4])
print('walls_LA_primary faces',len(W),'wall-only',len(wonly))
for k in ([idx[-1]] if len(sys.argv)<3 else [i for i in idx]):
    pass
P=X0+x.var(idx[-1],'displacement')
base={(a,b,j) for a,b,j,d in crossings(X0,wt,edges)}
cr=[c for c in crossings(P,wt,edges) if c[:3] not in base]
d=np.array([c[3] for c in cr]) if cr else np.zeros(0)
print('curve skin edges through wall-only faces at last state:',len(cr),' >0.5mm',int((d>0.5).sum()),' max %.2f'%(d.max() if len(d) else 0), ' at rest',len(base))
for a,b,j,dd in sorted(cr,key=lambda c:-c[3])[:6]: print('  %.2f mm edge rest'%dd, np.round((X0[a]+X0[b])/2,1))
