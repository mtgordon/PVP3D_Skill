"""Contact pressure / gap per plotted contact surface at a few converged states: which pairs are touching, and how
many facets. usage: py -3.10 contact_state.py RUN"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xplt import Xplt
run = sys.argv[1]
x = Xplt(os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs'), run, run + '.xplt'))
print('surface vars:', [it.get('name') for it in x.dict.get('surface', [])])
print('surfaces:', [s.get('name') if isinstance(s, dict) else s for s in getattr(x, 'surfaces', [])][:12])
conv = [i for i, s in enumerate(x.states) if s[1] == 0]
for k in [conv[int(len(conv)*f)] for f in (0.5, 0.8, 0.9)] + conv[-3:]:
    t = x.states[k][0]
    cp = x.var(k, 'contact pressure')
    gap = x.var(k, 'contact gap')
    if isinstance(cp, dict):
        out = []
        for sid in sorted(cp):
            a = np.asarray(cp[sid]).ravel(); g = np.asarray(gap[sid]).ravel() if isinstance(gap, dict) and sid in gap else np.array([0])
            out.append(f's{sid}: n>0 {int((a > 1e-9).sum())}/{a.size}, max p {a.max():.2e}, min gap {g.min():.2f}')
        print(f't {t:.4f}: ' + '; '.join(out))
    else:
        print(t, type(cp))
