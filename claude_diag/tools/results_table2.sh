#!/bin/bash
# Session 2026-09-23 results: markdown rows (run | base | change | t reached | retries | end | failure location)
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
LIVE=$(powershell -NoProfile -Command 'Get-CimInstance Win32_Process -Filter "Name='"'"'febio4.exe'"'"'" | ForEach-Object { $_.CommandLine }' 2>/dev/null)
row(){ n="$1"; b="$2"; d="$3"; w="$4"; L="$R/$n/$n.log"
  t=$(grep '^------- converged at time' "$L" 2>/dev/null | tail -1 | awk '{print $NF}')
  r=$(grep -c 'failed to converge' "$L" 2>/dev/null)
  s="stopped (by me)"; echo "$LIVE" | grep -q "$n.feb" && s="running"
  grep -q "N O R M A L   T E R M" "$L" 2>/dev/null && s="finished"; grep -q "E R R O R   T E R M" "$L" 2>/dev/null && s="failed"
  printf "| %s | %s | %s | %s | %s | %s | %s |\n" "$n" "$b" "$d" "${t:-0}" "${r:-0}" "$s" "$w"; }
echo "| Run | Base | Change (one per row) | t reached | retries | end | failure location |"
echo "|---|---|---|---|---|---|---|"
row K0_control TRUSS_fixed "none (control rerun)" "old wall; soft-mode stall (see G0)"
row P1_pvwla TRUSS_fixed "+ Abaqus Int-PVW-LA contact (sliding-elastic, offset 2 mm)" "identical to K0 (contact not engaged: gap 2.48 mm > 2 mm)"
row C1_pebavw TRUSS_fixed "+ 208 PeB-bottom facets in AVW-PVW contact (as Abaqus)" "identical to K0"
row X1_extconst TRUSS_fixed "connector springs extend=constant (Abaqus default)" "start-up retries (zero tangent), then old wall"
row L1_laload TRUSS_fixed "+ Abaqus Load-LA pressure (0.014 MPa, all LA)" "LA bulges 17-25 mm at 25% load; see LLd"
row F2_usll_snn0 TRUSS_fixed "USL-L_fan element normals (shell_normal_nodal=0)" "old wall"
row F1_fans_snn0 TRUSS_fixed "all 14 fans element normals" "old wall, slightly later"
row G0_diag TRUSS_fixed "diagnostic: 1 thread, negJ console, per-iteration plot" "divergence starts at L fan attachments: PVW node 2346 (USL-L tie surface), AVW-Para-L node 16287 (tnof); negJ AVW 1/1213/1303, USL-L_fan 12624"
row B1_vwmid TRUSS_fixed "+ Abaqus BC-VW-mid XSYMM (ux=0, 328 midline nodes)" "0 retries to 0.8587; then USL-L_fan 13118 / hex 7115 most strained"
row B1P1_vwmid_pvwla B1 "+ Abaqus Int-PVW-LA contact" "identical to B1 (0.858674); contact not engaged"
row B1C1_vwmid_pebavw B1 "+ 208 PeB-bottom facets" "identical to B1 (0.858674)"
row B1X1_vwmid_extconst B1 "connector springs extend=constant" "start-up retries, then B1 wall"
row B2_vwmid_nomidsprings B1 "- 328 midline stab springs (not in Abaqus)" "soft-mode stall at 0.504 (springs still hold y/z)"
row B1S7_diag_hex7115x10 B1 "DIAGNOSTIC: sliver hex 7115 10x stiffer" "B1 wall unchanged: hex 7115 is not the cause"
row B1F1_vwmid_fans_snn0 B1 "all 14 fans element normals" "passes B1 wall; crawls at 99.7% load (AVW-Para-R apex 23555, PM_PeB_Right 9358)"
row B1L1_vwmid_laload B1 "+ Abaqus Load-LA pressure" "early: LA_ICM node 23441 panel + USL-L_fan negJ (B1L1d)"
row B1L1d_diag B1L1 "diagnostic to t=0.1: 1 thread, negJ, per-iteration plot" "LA_ICM node 23441 (mid-panel), PVW 2347-2350 at P-arcus-L/USL-L; negJ USL-L_fan 12833/12803/13118"
row B1F1L1_fans_laload B1F1 "+ Abaqus Load-LA pressure" "LA-pressure wall (as L1)"
row LLd_diag_laload_wall B1F1L1 "diagnostic to t=0.35: 1 thread, negJ, per-iteration plot" "divergence starts at Posterior_Arcus_Right chain node 23723 (tied to LA 22542/23002/23008) and LA_PCMPRM node 22543; negJ then in AVW/PVW/PeB hexes"
row ALL1_faithful B1P1 "+ PeB facets + connectors constant" "B1 wall (no fan fix)"
row ALL2_faithful_fans ALL1 "+ all fans element normals  = TRUSS_fixed2" "crawls at 99.7% load, like B1F1"
row N1B1_noLA_vwmid "N1 (no LA)" "+ Abaqus BC-VW-mid XSYMM" "none: runs to full load"
