import os, sys, re, numpy as np
from collections import Counter, defaultdict
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load
run, out = sys.argv[1], sys.argv[2]
t, X, E = load(run); X = {k: np.asarray(v) for k,v in dict(X).items()}
HF8=[(0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]; HF6=[(0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)]
# PM_Plane boundary edges
pq=[c for n,c in E.values() if n=='PM_Plane']
ec=Counter(frozenset((q[i],q[(i+1)%4])) for q in pq for i in range(4))
bnd=[tuple(e) for e,n in ec.items() if n==1]
def on_edge(p):
    best=None
    for a,b in bnd:
        A,B=X[a],X[b]; s=np.clip(np.dot(p-A,B-A)/np.dot(B-A,B-A),0,1); q=A+s*(B-A); d=np.linalg.norm(q-p)
        if best is None or d<best[0]: best=(d,q,(a,b))
    return best
an=[int(v) for v in re.findall(r'\d+', re.search(r'<NodeSet name="curve_end_PM_anchors">([^<]*)</NodeSet>', t).group(1))]
print('anchor distance to PM_Plane boundary:', max(on_edge(X[a])[0] for a in an).round(2))
props=[]
for sgn,corner,avw,anc in ((-1,25896,1775,27597),(1,None,None,27600)):
    if corner is None:
        # mirror on the right: find strap/curve node near mirrored left
        def near(p, pool): return min(pool, key=lambda v: np.linalg.norm(X[v]-p))
        sr={v for n,c in E.values() if n=='PM_avw_bottom_right_sheet_solid' for v in c}
        cv={v for n,c in E.values() if n.startswith('canal_wrap') for v in c}
        m=lambda v: X[v]*np.array([-1,1,1])
        corner=near(m(25896)+np.array([2.5,0,0]), sr&cv); avw=near(m(1775)+np.array([2.5,0,0]), cv)
    V=X[anc]-X[corner]; tgt=X[avw]+V; d,q,e=on_edge(tgt)
    print('side',sgn,'corner',corner,X[corner].round(2),'avw',avw,X[avw].round(2),'anchor',anc,X[anc].round(2),'arch target',tgt.round(2),'-> on PM_Plane edge',q.round(2),'off',round(d,2), 'curve edge +',round(np.linalg.norm(X[avw]-X[corner]),2),'arch edge +',round(np.linalg.norm(q-X[anc]),2))
    props.append([X[corner],X[avw],q,X[anc]])
fig,axs=plt.subplots(1,2,figsize=(15,7))
for ax,(xl,yl) in zip(axs,(((-45,45),(-5,40)),((-24,-1),(8,32)))):
    for k,(n,c) in E.items():
        P=np.array([X[v] for v in c]); 
        if P[:,2].mean()<-60 or P[:,2].mean()>-25 or n.startswith('LA'): continue
        col='orange' if 'sheet_solid' in n else 'mediumturquoise' if n.startswith('canal_wrap') else 'yellowgreen' if n=='PM_Plane' else 'olive' if n=='_PickedSet347' else None
        if col is None: continue
        fl=HF8 if len(c)==8 else HF6 if len(c)==6 else [tuple(range(len(c)))]
        for f in fl[:2] if len(c)==8 else fl[:1]:
            ax.add_patch(plt.Polygon(P[list(f)][:,:2],fc=col,ec='k',lw=0.1,alpha=0.8))
    for p in props:
        ax.add_patch(plt.Polygon(np.array(p)[:,:2],fc='red',ec='darkred',lw=1,alpha=0.55))
    ax.plot([X[a][0] for a in an],[X[a][1] for a in an],'k.',ms=3)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.invert_yaxis(); ax.set_aspect('equal'); ax.grid(lw=0.2)
axs[0].set_title('proposal (red): strap widened to the curve element at the AVW + more PM arch\nseen from distal; orange strap, turquoise curve, olive AVW, green PM_Plane, black = fixed strap anchors',fontsize=9)
axs[1].set_title('left side, zoom',fontsize=9)
plt.tight_layout(); plt.savefig(out,dpi=120)
