"""curve_strap primary nodes near / behind the strap faces at rest (2026-10-22).
usage: py -3.10 strap_rest_gap.py RUN"""
import os, sys, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load
t,X,E=load(sys.argv[1]); X={k:np.asarray(v) for k,v in dict(X).items()}
sec=[list(map(int,f.split(','))) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="curve_strap_secondary">(.*?)</Surface>',t,re.S).group(1))]
prim={int(v) for f in re.findall(r'>([\d,]+)</',re.search(r'<Surface name="curve_strap_primary">(.*?)</Surface>',t,re.S).group(1)) for v in f.split(',')}
res=[]
for v in prim:
    p=X[v]; best=9
    for f in sec:
        if v in f: continue
        A,B,C=[X[q] for q in f[:3]]; n=np.cross(B-A,C-A); n/=np.linalg.norm(n); d=np.dot(p-A,n)
        q=p-d*n; v0,v1,v2=B-A,C-A,q-A; d00,d01,d11,d20,d21=v0@v0,v0@v1,v1@v1,v2@v0,v2@v1; den=d00*d11-d01*d01
        l1=(d11*d20-d01*d21)/den; l2=(d00*d21-d01*d20)/den
        if l1>=-0.05 and l2>=-0.05 and l1+l2<=1.05 and abs(d)<abs(best): best=d
    if best<9: res.append((best,v))
res.sort(key=lambda r: abs(r[0]))
print('curve_strap primary nodes within 0.6 mm of a strap face at rest (signed distance, + = outside):', [(round(d,3),v) for d,v in res if abs(d)<0.6][:20])
