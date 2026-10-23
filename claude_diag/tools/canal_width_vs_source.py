"""Canal width (walls + curves + cervix) vs the source AVW/PVW/cervix by level (2026-10-22).
usage: py -3.10 canal_width_vs_source.py RUN"""
import os, sys, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load
L=open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'PVP3DModel_job.inp')).read().split('\n')[3536:19640]
X={}; El={}; sets={}; i=0; mode=None; cur=None; gen=False
for ln in L:
    if ln.startswith('*'):
        l=ln.lower(); mode=None
        if l.startswith('*node'): mode='n'
        elif l.startswith('*element'): mode='e'
        elif l.startswith('*elset'):
            mode='s'; cur=re.search(r'elset=([^,\s]+)',ln).group(1); gen='generate' in l; sets[cur]=[]
        continue
    if not ln.strip(): continue
    v=[s.strip() for s in ln.split(',') if s.strip()]
    if mode=='n': X[int(v[0])]=np.array(list(map(float,v[1:4])))
    elif mode=='e': El[int(v[0])]=list(map(int,v[1:9]))
    elif mode=='s':
        if gen: a,b,s_=map(int,v[:3]); sets[cur]+=list(range(a,b+1,s_))
        else: sets[cur]+=list(map(int,v))
def nodes(es): return np.array([X[n] for e in es for n in El[e]])
src={'AVW':nodes(sets['_PickedSet347']),'PVW':nodes(sets['_PickedSet64']),'cervix':nodes(sets['_PickedSet346'])}; SALL=np.vstack(list(src.values()))
t,FX,E=load(sys.argv[1]); FX={k:np.asarray(v) for k,v in dict(FX).items()}
def fn(pred): return np.array([FX[v] for k,(n,c) in E.items() if pred(n) for v in c])
our={'AVW':fn(lambda n:n=='_PickedSet347'),'PVW':fn(lambda n:n=='_PickedSet64'),'curves':fn(lambda n:n.startswith('canal_wrap_hex')),'cervix':fn(lambda n:n=='_PickedSet346')}
allour=np.vstack([our['AVW'],our['PVW'],our['curves']])
print('z band | source AVW x | source PVW x | src width | ours AVW x | ours PVW x | ours curves x | our width (walls+curves) | ratio')
for z0 in range(-42,36,4):
    def rng(A):
        s=A[(A[:,2]>=z0)&(A[:,2]<z0+4)]
        return (s[:,0].min(),s[:,0].max()) if len(s) else None
    r={k:rng(v) for k,v in list(src.items())}; o={k:rng(v) for k,v in our.items()}
    sw=[r[k] for k in ('AVW','PVW') if r[k]]; ow=[o[k] for k in ('AVW','PVW','curves') if o[k]]
    if not sw or not ow: continue
    swd=max(b for a,b in sw)-min(a for a,b in sw); owd=max(b for a,b in ow)-min(a for a,b in ow)
    f=lambda p: f'{p[0]:6.1f}..{p[1]:5.1f}' if p else '      -      '
    print(f'{z0:4d}..{z0+4:<4d}| {f(r["AVW"])} | {f(r["PVW"])} | {swd:5.1f} | {f(o["AVW"])} | {f(o["PVW"])} | {f(o["curves"])} | {owd:5.1f} | {owd/swd:.2f}')
print('source AVW y range', src['AVW'][:,1].min().round(1), src['AVW'][:,1].max().round(1), ' ours', our['AVW'][:,1].min().round(1), our['AVW'][:,1].max().round(1))
print('--- with the cervix: source AVW+PVW+cervix vs ours AVW+PVW+cervix+curves')
OALL=np.vstack([our['AVW'],our['PVW'],our['curves'],our['cervix']])
for z0 in range(-42,46,4):
    s=SALL[(SALL[:,2]>=z0)&(SALL[:,2]<z0+4)]; o=OALL[(OALL[:,2]>=z0)&(OALL[:,2]<z0+4)]
    c=src['cervix'][(src['cervix'][:,2]>=z0)&(src['cervix'][:,2]<z0+4)]; oc=our['cervix'][(our['cervix'][:,2]>=z0)&(our['cervix'][:,2]<z0+4)]
    if len(s) and len(o): print(f'{z0:4d}..{z0+4:<4d} src {s[:,0].min():6.1f}..{s[:,0].max():5.1f} ({np.ptp(s[:,0]):5.1f})  ours {o[:,0].min():6.1f}..{o[:,0].max():5.1f} ({np.ptp(o[:,0]):5.1f})  ratio {np.ptp(o[:,0])/np.ptp(s[:,0]):.2f}   cervix src {np.ptp(c[:,0]) if len(c) else 0:5.1f} ours {np.ptp(oc[:,0]) if len(oc) else 0:5.1f}')
