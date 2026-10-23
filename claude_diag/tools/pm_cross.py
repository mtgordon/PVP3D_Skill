"""Canal surface edges (curves / AVW / PVW / PeB) through the PM_Plane arch (2026-10-22/23).
usage: py -3.10 pm_cross.py RUN T1,T2,..."""
import sys, os, numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from curve_cross import crossings, RUNS, FI
from strap_cross import clusters
from xplt_reader import Xplt, converged_state_indices
run=sys.argv[1]
x=Xplt(os.path.join(RUNS,run,run+'.xplt')); idx=converged_state_indices(x,os.path.join(RUNS,run,run+'.log'))
print([ (d.get('name'), d['etype'] if 'etype' in d else '', d['conn'].shape) for d in x.domains if (d.get('name') or '') in ('PM_Plane',)])
pm=[d for d in x.domains if d['conn'].shape==(70,4)][0]['conn'][:,:4]
pt=np.array([q[[0,1,2]] for q in pm]+[q[[0,2,3]] for q in pm])
def ext(prefixes):
    H=np.vstack([d['conn'][:,:8] for d in x.domains if any((d.get('name') or '').startswith(p) for p in prefixes) and d['conn'].shape[1]>=8])
    cnt=Counter(tuple(sorted(h[list(f)])) for h in H for f in FI)
    quads=[h[list(f)] for h in H for f in FI if cnt[tuple(sorted(h[list(f)]))]==1]
    return np.array(sorted({tuple(sorted((q[k],q[(k+1)%4]))) for q in quads for k in range(4)}))
X0=x.X; tt=np.array([x.states[i][0] for i in idx])
for lab,pre in (('curves',('canal_wrap',)),('AVW',('_PickedSet347',)),('PVW',('_PickedSet64',)),('PeB',('_PickedSet66',))):
    ed=ext(pre); base={c[:3] for c in crossings(X0,pt,ed)}
    for T in sys.argv[2].split(','):
        i=idx[int(np.argmin(abs(tt-float(T))))]; P=X0+x.var(i,'displacement')
        cr=[c for c in crossings(P,pt,ed) if c[:3] not in base]
        d=[c[3] for c in cr]
        print(f'{lab:6s} t {x.states[i][0]:.3f}: edges through PM_Plane {len(cr)}, >0.5 mm {sum(v>0.5 for v in d)}, max {max(d) if d else 0:.2f}; rest-crossings {len(base)}; clusters', '; '.join(f'{n} at {tuple(cc)} {dm:.2f}' for n,cc,dm in clusters(X0,cr)[:4]))
