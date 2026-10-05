# usage: py -3.10 domain_mass.py MODEL.feb   (mass per domain: density x volume, grams; mass-only BeamDomains show 0 here)
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import xml.etree.ElementTree as ET
from febmodel import Feb
f = Feb(sys.argv[1])
root = ET.parse(sys.argv[1]).getroot()
mats = {m.get('name'): m for m in root.find('Material')}
md = root.find('MeshDomains')
X = f.nodes
out = []
for d in md:
    name, mat = d.get('name'), d.get('mat')
    m = mats.get(mat)
    rho = float(m.find('density').text) if m is not None and m.find('density') is not None else 0.0
    typ, conn = f.elem_blocks.get(name, (None, {}))
    vol = 0.0
    if d.tag == 'ShellDomain':
        t = None
        sh = d.find('shell_thickness')
        t = float(sh.text.split(',')[0]) if sh is not None else None
        # area of tri/quad elements
        for e, c in conn.items():
            P = np.array([X[n] for n in c])
            if len(c) == 3:
                a = 0.5*np.linalg.norm(np.cross(P[1]-P[0], P[2]-P[0]))
            else:
                a = 0.5*np.linalg.norm(np.cross(P[2]-P[0], P[3]-P[1]))
            vol += a*(t or 0)
    elif d.tag == 'SolidDomain':
        for e, c in conn.items():
            P = np.array([X[n] for n in c])
            if len(c) == 8:
                # hex volume by 5-tet decomposition approx via bounding: use 6 tets
                tets = [(0,1,3,4),(1,2,3,6),(1,4,5,6),(3,4,6,7),(1,3,4,6)]
                for a,b,cc,dd in tets:
                    vol += abs(np.dot(P[b]-P[a], np.cross(P[cc]-P[a], P[dd]-P[a])))/6
            elif len(c) == 4:
                vol += abs(np.dot(P[1]-P[0], np.cross(P[2]-P[0], P[3]-P[0])))/6
    out.append((name, d.tag, mat, rho, vol, rho*vol*1e6))  # grams (t -> g: 1e6)
for r in sorted(out, key=lambda r: -r[5]):
    print(f'{r[0]:32s} {r[1]:12s} {r[2][:30]:30s} rho {r[3]:.3g} vol {r[4]:9.1f} mm3 mass {r[5]:9.2f} g')
