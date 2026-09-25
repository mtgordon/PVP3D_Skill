"""Map Abaqus part-local nodes to FEBio node ids by coordinate."""
import re
import numpy as np
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import INP_FILE

INP = INP_FILE


def part_nodes(part, inp=INP):
    """Return {local_id: xyz} for *Part name=part."""
    lines = open(inp, encoding='latin-1').read().splitlines()
    out = {}
    inpart = False
    innode = False
    for ln in lines:
        s = ln.strip()
        if s.lower().startswith('*part,'):
            m = re.search(r'name=("?)(.+?)\1\s*$', s)
            inpart = m and m.group(2) == part
            innode = False
            continue
        if not inpart:
            continue
        if s.lower().startswith('*end part'):
            break
        if s.startswith('*'):
            innode = s.lower().startswith('*node') and not s.lower().startswith('*node output')
            continue
        if innode and s:
            v = [t.strip() for t in s.split(',')]
            out[int(v[0])] = np.array([float(v[1]), float(v[2]), float(v[3])])
    return out


def nearest_feb(feb, xyz, candidates=None):
    ids = np.array(sorted(candidates if candidates is not None else feb.nodes.keys()))
    P = np.array([feb.nodes[i] for i in ids])
    d = np.linalg.norm(P - xyz, axis=1)
    k = int(np.argmin(d))
    return int(ids[k]), float(d[k])
