"""Batch 16f: the stiff arcus beam (16d) broke the start-up at the LA-PeB sphincter connectors. One change each
from LP2: a 10x softer beam (G*As 2 N > the ~1.2 N compression, EI 4 N mm^2), or the same beam with pinned ends."""
import os
from variants4 import Model4, JOBS, emit

LP2 = os.path.join(JOBS, 'claude_diag', 'runs', 'LP2_yeoh_k1', 'LP2_yeoh_k1.feb')
m = Model4(LP2); m.add_arcus_beam(E=1.0, G=0.4); emit('LPBs_beam_soft', m)
m = Model4(LP2); m.add_arcus_beam(clamp_ends=False); emit('LPBp_beam_pinned', m)
