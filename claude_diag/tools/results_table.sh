#!/bin/bash
# name|description ; prints markdown rows with last converged t, retries, status
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
J="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LIVE=$(powershell -NoProfile -Command 'Get-CimInstance Win32_Process -Filter "Name='"'"'febio4.exe'"'"'" | ForEach-Object { $_.CommandLine }' 2>/dev/null)
row(){ L="$1"; d="$2"; n="$3"
  t=$(grep '^------- converged at time' "$L" 2>/dev/null | tail -1 | awk '{print $NF}')
  r=$(grep -c 'failed to converge' "$L" 2>/dev/null)
  s="stopped (by me)"; echo "$LIVE" | grep -q "$(basename "$L" .log).feb" && s="running"; grep -q "N O R M A L   T E R M" "$L" && s="finished"; grep -q "E R R O R   T E R M" "$L" && s="failed"
  printf "| %s | %s | %s | %s | %s |\n" "$n" "$d" "${t:-0}" "$r" "$s"; }
echo "| Run | Change (single step from the row it builds on) | t reached | retries | end |"
echo "|---|---|---|---|---|"
row $J/run_truss/PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.log "your TRUSS file" "TRUSS (yours)"
row $J/run_truss_notie/PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_NOTIE.log "your TRUSS, ties off" "TRUSS_NOTIE (yours)"
row $J/jobs/PVP3DModel_v21_v3_HYBRID_withLA_CORRECTED_testing.log "your best LA run (ties inert)" "CORRECTED_testing (yours)"
row $R/E3_v32_corrctrl/E3_v32_corrctrl.log "no-LA v32 + LA solver settings" "E3"
row $R/E0_base/E0_base.log "CORRECTED, unchanged (control)" "E0"
row $R/E4_base_v32ctrl/E4_base_v32ctrl.log "E0 with v32 solver settings" "E4"
row $R/E1_unmerge/E1_unmerge.log "E0, LA un-merged from sphincter fans (LA unloaded)" "E1"
row $R/E2_abqconn/E2_abqconn.log "E1 + 8 Abaqus side connectors" "E2"
row $R/E5_tubes_connLRP/E5_tubes_connLRP.log "E1 + all 40 Abaqus connectors" "E5"
row $R/T1_truss_lc/T1_truss_lc.log "TRUSS, ties -> exact linear constraints (fan merge kept)" "T1"
row $R/T2_truss_lc_connLR/T2_truss_lc_connLR.log "T1, un-merged + 16 side connectors" "T2"
row $R/T3_truss_lc_connLRP/T3_truss_lc_connLRP.log "T2 + 24 posterior connectors" "T3"
row $R/T3s_truss_lc_connLRP_pins/T3s_truss_lc_connLRP_pins.log "T3 + remesh-lost pins restored" "T3s"
row $R/T3f_faithful/T3f_faithful.log "T3s - redundant sphincter fans" "T3f"
row $R/T3f_k/T3f_k.log "T3f + LA k=0.349" "T3f_k"
row $R/T3f_3f/T3f_3f.log "T3f + three-field LA shells" "T3f_3f"
row $R/T3f_t15/T3f_t15.log "T3f, LA t=1.5 membrane-equivalent" "T3f_t15"
row $R/D1_faithful_dyn/D1_faithful_dyn.log "T3f, implicit DYNAMIC + Abaqus densities" "D1"
row $R/T3f_snn0/T3f_snn0.log "**T3f + LA element normals (shell_normal_nodal=0)**" "**T3f_snn0**"
row $R/T3f_k_snn0/T3f_k_snn0.log "**T3f_snn0 + LA k=0.349 -> TRUSS_fixed**" "**T3f_k_snn0**"
row $R/E0n_base_snn0/E0n_base_snn0.log "E0 + LA element normals only (fan merge kept)" "E0n"
row $R/N1_v32_pins/N1_v32_pins.log "no-LA v32 + remesh-lost pins restored" "N1"
