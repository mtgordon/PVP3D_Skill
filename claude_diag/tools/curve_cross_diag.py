"""Which surfaces take part in the deep (> 0.5 mm) curve crossings at the last state, and whether they are in the curve self-contact (curve_all_self): "PS" = in both primary and secondary, "-" = not in it. Rest gap = distance between the edge and the face at rest. (2026-10-25, the colleague.)
usage: py -3.10 curve_cross_diag.py RUN"""
import os,sys,numpy as np,xml.etree.ElementTree as ET
sys.path.insert(0,'/home/claude/pvp3d_skill/claude_diag/tools')
from curve_cross import surface, crossings, RUNS
from paths import SKILL_SCRIPTS; sys.path.insert(0,SKILL_SCRIPTS)
from xplt_reader import Xplt, converged_state_indices
run=sys.argv[1]
x=Xplt(os.path.join(RUNS,run,run+'.xplt')); idx=converged_state_indices(x,os.path.join(RUNS,run,run+'.log'))
tris,edges=surface(x); X0=x.X
# xplt node index -> feb node id: assume order of ids sorted
r=ET.parse(os.path.join(RUNS,run,run+'.feb')).getroot(); m=r.find('Mesh')
ids=[int(n.get('id')) for b in m.findall('Nodes') for n in b]
surf={s.get('name'):{int(v) for f in s for v in f.text.split(',')} for s in m.findall('Surface')}
pair=[p for p in m.findall('SurfacePair') if p.get('name')=='curve_all_self'][0]
prim=surf[pair.find('primary').text]; sec=surf[pair.find('secondary').text]
H=np.vstack([d['conn'][:,:8] for d in x.domains if (d.get('name') or '').startswith('canal_wrap_hex')])
node2el={}
for k,h in enumerate(H):
    for n in h: node2el.setdefault(int(n),set()).add(k)
P=X0+x.var(idx[-1],'displacement')
cr=crossings(P,tris,edges)
base={(a,b,j) for a,b,j,d in crossings(X0,tris,edges)}
cr=[c for c in cr if c[:3] not in base and c[3]>0.5]
print(run,'deep crossings (>0.5 mm):',len(cr))
for a,b,j,d in sorted(cr,key=lambda c:-c[3])[:30]:
    t=tris[j]; e_el=node2el[int(a)]|node2el[int(b)]; t_el=set.intersection(*[node2el[int(v)] for v in t]) if True else set()
    t_el2=set().union(*[node2el[int(v)] for v in t])
    share = 'same/adjacent hex' if e_el & t_el2 else 'separate'
    rest=np.linalg.norm((X0[a]+X0[b])/2-X0[t].mean(0))
    fid=lambda n: ids[int(n)]
    inS=lambda ns: ('P' if all(fid(n) in prim for n in ns) else '')+('S' if all(fid(n) in sec for n in ns) else '')
    print('%.2f mm  edge rest (%.1f,%.1f,%.1f) %s | tri %s | rest gap %.1f mm | contact sets edge:%s tri:%s'%(d,*((X0[a]+X0[b])/2),share,'',rest,inS([a,b]) or '-',inS(t) or '-'))
