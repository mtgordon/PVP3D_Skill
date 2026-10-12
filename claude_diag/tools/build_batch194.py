"""Batch 194 (2026-10-06, the user: "yes, start the distal softening test"). The tube's distal end does not peel back
(evert) as much as the user expects. A LOCAL soft zone, not whole tissue: on L187_selfall_botfull_P1 (the tube with the
user's equal PM-arch springs, no curve pressure, ran to t 1), the AVW, PVW and curve elements whose centroid lies within
ZONE mm of the hymen at rest (height > -ZONE; 0 = hymen, PM_Plane fit as pop_q_bump.py) get their own domains
<name>_dist with copies of their materials at c1, c2, k x FACTOR (NOT IN SOURCE). Their walls are x0.35 already, so the
zone is at x0.175 of the source. The new domains keep the walls' mass damping (C 200).
  L194_botfull_distsoft05_P1   ZONE 15 mm, FACTOR 0.5
The POP-Q tools read <name>_dist as <name>.
usage: py -3.10 build_batch194.py"""
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
SRC = 'L187_selfall_botfull_P1'
NAME = 'L194_botfull_distsoft05_P1'
ZONE, FACTOR = 15.0, 0.5
DOMS = {'_PickedSet347': 'Vagina_AVW', '_PickedSet64': 'Vagina_PVW', 'canal_wrap_hex_avw': 'Vagina_AVW',
        'canal_wrap_hex_pvw': 'Vagina_PVW'}

if __name__ == '__main__':
    t = open(os.path.join(RUNS, SRC, SRC + '.feb'), encoding='utf-8').read()
    mesh = ET.fromstring(t).find('Mesh')
    X = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in mesh.findall('Nodes') for n in b}
    ref = ET.parse(os.path.join(RUNS, 'L26_springs_la3_pm', 'L26_springs_la3_pm.feb')).getroot().find('Mesh')
    XR = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for b in ref.findall('Nodes') for n in b}
    P = np.array([XR[i] for i in range(23817, 24213)])
    c0 = P.mean(0)
    nrm = np.linalg.svd(P - c0)[2][-1]
    cx = [int(v) for e in mesh.findall('Elements') if e.get('name') == '_PickedSet346' for el in e for v in el.text.split(',')]
    nn = nrm if np.dot(np.mean([X[v] for v in cx], 0) - c0, nrm) < 0 else -nrm
    log, newblocks, newdoms, newdamp = [], '', '', ''
    for dom, mat in DOMS.items():
        a = t.index(f'<Elements type="hex8" name="{dom}">')
        b = t.index('</Elements>', a)
        lines = t[a:b].split('\n')
        head, body = lines[0], [l for l in lines[1:] if '<elem' in l]
        keep, move = [], []
        for l in body:
            ids = [int(v) for v in re.search(r'>([^<]+)<', l).group(1).split(',')]
            (move if (np.mean([X[i] for i in ids], 0) - c0) @ nn > -ZONE else keep).append(l)
        t = t[:a] + '\n'.join([head] + keep) + '\n\t\t' + t[b:]
        newblocks += f'\t\t<Elements type="hex8" name="{dom}_dist">\n' + '\n'.join(move) + '\n\t\t</Elements>\n'
        newdoms += f'\t\t<SolidDomain name="{dom}_dist" mat="{mat}_dist" />\n'
        newdamp += (f'\t\t<body_load name="wall_damping_{dom.strip("_")}_dist" type="mass damping" elem_set="{dom}_dist">\n'
                    f'\t\t\t<C>200</C>\n\t\t</body_load>\n')
        log.append(f'{dom}: {len(move)} of {len(body)} elements -> {dom}_dist')
    i = t.index('\t</Mesh>')
    t = t[:i] + newblocks + t[i:]
    i = t.index('\t</MeshDomains>')
    t = t[:i] + newdoms + t[i:]
    i = t.index('\t</Loads>')
    t = t[:i] + newdamp + t[i:]
    mats = ''
    nid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t))
    for mat in ('Vagina_AVW', 'Vagina_PVW'):
        a = t.index(f'name="{mat}" type="Yeoh">')
        a = t.rindex('<material', 0, a)
        b = t.index('</material>', a) + len('</material>')
        blk = t[a:b]
        for k in ('c1', 'c2', 'k'):
            blk = re.sub(rf'<{k}>([^<]+)</{k}>', lambda m: f'<{k}>{float(m.group(1)) * FACTOR:.7g}</{k}>', blk)
        nid += 1
        blk = re.sub(r'id="\d+"', f'id="{nid}"', blk, count=1).replace(f'name="{mat}"', f'name="{mat}_dist"')
        mats += '\t\t' + blk + '\n'
    i = t.index('\t</Material>')
    t = t[:i] + mats + t[i:]
    write(NAME, t, f'base: claude_diag/runs/{SRC}/{SRC}.feb (tools/build_batch194.py)\n'
                   f'NOT IN SOURCE (a LOCAL soft zone, not whole tissue): the AVW, PVW and curve elements within {ZONE:g} mm of the '
                   f'hymen at rest moved to <name>_dist domains with Vagina_AVW_dist / Vagina_PVW_dist = c1, c2, k x {FACTOR:g} '
                   f'(walls already x0.35, so x{0.35 * FACTOR:g} of the source); damping C 200 kept. ' + '; '.join(log) + '\n')
    print('\n'.join(log))
