"""Every Abaqus CONN3D2 connector: end parts, behavior family, count. Ends given as instance.node or
assembly-level node ids (reference points)."""
import re
from collections import Counter, defaultdict
from abq_surf import lines

L = lines()
rows = []
for i, ln in enumerate(L):
    if ln.strip().lower().startswith('*element, type=conn3d2'):
        data = L[i + 1].strip()
        sec = L[i + 2].strip()
        v = [t.strip() for t in data.split(',')]
        a, b = v[1], v[2] if len(v) > 2 else '(ground)'
        m = re.search(r'behavior=([^,\s]+)', sec)
        rows.append((a, b, m.group(1) if m else '?'))
def part(e):
    return e.split('.')[0] if '.' in e else 'ASSEMBLY-node'
fam = lambda s: re.sub(r'-(I|II|III|IV)-.*', '', s).replace('-x%stiff', '')
c = Counter((part(a), part(b), fam(s)) for a, b, s in rows)
print(f'{len(rows)} CONN3D2 connectors')
for (pa, pb, f), n in sorted(c.items(), key=lambda kv: (-kv[1])):
    print(f'  {n:4d}  {pa:22s} -> {pb:22s}  {f}')
