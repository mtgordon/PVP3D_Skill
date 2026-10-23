"""Largest node move between consecutive plot states over a run's last N states (default 30): spots spurious jumps before
an end state is used. Prints per step: time, max |du| (mm), the node (1-based plot index) and the 95th percentile.
usage: py -3.10 last_moves.py RUN [N]"""
import os
import sys

import numpy as np

from paths import SKILL_SCRIPTS  # noqa: E402
sys.path.insert(0, SKILL_SCRIPTS)
from xplt_reader import Xplt  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')
run = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
x = Xplt(os.path.join(RUNS, run, run + '.xplt'))
ids = list(range(max(0, len(x.states) - n - 1), len(x.states)))
prev = None
worst = (0, None, None)
for i in ids:
    u = np.asarray(x.var(i, 'displacement'), dtype=float)
    if prev is not None:
        d = np.linalg.norm(u - prev, axis=1)
        j = int(d.argmax())
        print(f't {x.states[i][0]:.4f}  max {d[j]:7.3f} mm (plot index {j + 1} = node id {int(x.node_ids[j])})  p95 {np.percentile(d, 95):.3f}')
        if d[j] > worst[0]:
            worst = (d[j], x.states[i][0], f'plot index {j + 1} = node id {int(x.node_ids[j])}')
    prev = u
print(f'{run}: largest move {worst[0]:.3f} mm at t {worst[1]} ({worst[2]}) over the last {len(ids) - 1} steps')
