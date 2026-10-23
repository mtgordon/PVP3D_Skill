"""Turn FEBio's per-element negative-jacobian console output into a domain-level report.

Produce the console output (FEBio 4; the plain `-i` run only prints a count):

    printf 'set output_negative_jacobians 1\\nrun -i model.feb\\nquit\\n' | \\
        OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 febio4.exe > model.console.txt

(single-threaded so the failure you diagnose is the same one every time), then:

    py -3 negjac_report.py model.feb model.console.txt

Prints each reported element with its domain, how often it was reported, the most
negative jacobian seen, and its node IDs. Elements that recur across retries and
across different variants are the ones to inspect (reference geometry first:
geometry-cross-referencing.md Rules 4 and 11).
Needs Python 3 + numpy; imports feb_model.py from this directory.
"""
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feb_model import Feb  # noqa: E402

feb_path, console_path = sys.argv[1], sys.argv[2]
txt = open(console_path, encoding='latin-1').read()
hits = re.findall(r'Negative jacobian was detected at element (\d+) at gauss point (\d+)\s*\*?\s*\n'
                  r'\s*\*?\s*jacobian = ([-0-9.eE+]+)', txt)
if not hits:
    sys.exit('no per-element reports found (was output_negative_jacobians set before run?)')
f = Feb(feb_path)
count, worst = defaultdict(int), {}
for eid, gp, jac in hits:
    e = int(eid)
    count[e] += 1
    worst[e] = min(worst.get(e, 0.0), float(jac))
by_dom = defaultdict(int)
for e, c in count.items():
    by_dom[f.elem_owner.get(e, '?')] += c
print('reports per domain:', dict(sorted(by_dom.items(), key=lambda kv: -kv[1])))
print(f'{"element":>8}  {"domain":28s} {"reports":>7}  {"worst J":>10}  nodes')
for e, c in sorted(count.items(), key=lambda kv: (-kv[1], worst[kv[0]]))[:30]:
    blk = f.elem_owner.get(e)
    conn = f.elem_blocks[blk][1][e] if blk else []
    print(f'{e:8d}  {str(blk):28s} {c:7d}  {worst[e]:10.4g}  {conn}')
