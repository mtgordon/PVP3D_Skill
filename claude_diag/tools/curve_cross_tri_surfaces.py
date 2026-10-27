"""Which feb surfaces the faces pierced in the deep (> 0.5 mm) curve crossings belong to (tied nodes, curve_all, curve_endfaces of a C018/C022 feb).
usage: py -3.10 curve_cross_tri_surfaces.py RUN FEB_WITH_curve_endfaces (2026-10-26, the colleague)"""
import os,sys,re,numpy as np,xml.etree.ElementTree as ET
sys.path.insert(0,'/home/claude/pvp3d_skill/claude_diag/tools')
from curve_cross import surface, crossings, RUNS
from paths import SKILL_SCRIPTS; sys.path.insert(0,SKILL_SCRIPTS)
from xplt_reader import Xplt, converged_state_indices
run=sys.argv[1]; feb2=sys.argv[2]
x=Xplt(os.path.join(RUNS,run,run+'.xplt')); idx=converged_state_indices(x,os.path.join(RUNS,run,run+'.log'))
tris,edges=surface(x); X0=x.X
t=open(os.path.join(RUNS,run,run+'.feb')).read()
ids=[int(v) for v in re.findall(r'<node id="(\d+)"',t)]
tied={int(v) for v in re.findall(r'<node id="(\d+)"',re.search(r'name="canal_wrap_ends".*?</constraint>',t,re.S).group(0))}
ca={int(v) for f in re.findall(r'>([^<]+)</quad4>',re.search(r'<Surface name="curve_all">(.*?)</Surface>',t,re.S).group(1)) for v in f.split(',')}
t2=open(feb2).read()
ef={int(v) for f in re.findall(r'>([^<]+)</quad4>',re.search(r'<Surface name="curve_endfaces">(.*?)</Surface>',t2,re.S).group(1)) for v in f.split(',')}
# which other surfaces contain these nodes
surfs={m.group(1):{int(v) for f in re.findall(r'>([^<]+)</(?:quad4|tri3)>',m.group(2)) for v in f.split(',')} for m in re.finditer(r'<Surface name="([^"]+)">(.*?)</Surface>',t,re.S)}
P=X0+x.var(idx[-1],'displacement')
base={(a,b,j) for a,b,j,d in crossings(X0,tris,edges)}
cr=[c for c in crossings(P,tris,edges) if c[:3] not in base and c[3]>0.5]
for a,b,j,d in sorted(cr,key=lambda c:-c[3])[:8]:
    tn=[ids[int(v)] for v in tris[j]]
    print('%.2f'%d,'tri nodes',tn,'tied',sum(n in tied for n in tn),'in curve_all',sum(n in ca for n in tn),'in endfaces',sum(n in ef for n in tn),
          'surfaces with all nodes:',[s for s,ns in surfs.items() if all(n in ns for n in tn)])
