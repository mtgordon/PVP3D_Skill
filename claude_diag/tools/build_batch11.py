"""Batch 11: fan shells with element normals (the USL-L_fan pinch found with shell_pinch_top.py)."""
from variants2 import Model2, BEST, emit, set_shell_normal_nodal, FANS

m = Model2(BEST); set_shell_normal_nodal(m, ('USL-L_fan',)); emit('F2_usll_snn0', m)
m = Model2(BEST); set_shell_normal_nodal(m, FANS); emit('F1_fans_snn0', m)
