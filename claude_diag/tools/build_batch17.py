"""Batch 17: the user's fan convention done visibly and softly. One change per step:
LPFm  = LP2 + P-arcus fans MERGED onto the arcus chain nodes (shared nodes; should reproduce LPF_fanties, the tie version)
LPFms = LPFm + soft fan material (Ogden fitted to the 13 connectors per side)
LPFms_DM1 = LP2DM_dyn_chainmass + merged soft fans (like-for-like with Abaqus; compare LPC_DM1_pconn, LPF_DM1_fanties);
            chain-mass trusses written as line2 (format only, so FEBio Studio can open it)."""
import os
from variants4 import Model4, JOBS, emit

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
LP2 = os.path.join(RUNS, 'LP2_yeoh_k1', 'LP2_yeoh_k1.feb')
DM1 = os.path.join(RUNS, 'LP2DM_dyn_chainmass', 'LP2DM_dyn_chainmass.feb')

m = Model4(LP2); m.merge_parcus_fans(); emit('LPFm_merged', m)
m = Model4(LP2); m.merge_parcus_fans(); m.set_parcus_fan_material(); emit('LPFms_merged_soft', m)
m = Model4(DM1); m.truss2_to_line2(); m.merge_parcus_fans(); m.set_parcus_fan_material(); emit('LPFms_DM1', m)
