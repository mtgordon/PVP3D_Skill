#!/bin/bash
# Summarize progress of every run under claude_diag/runs (or given names)
RUNS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
cd "$RUNS"
names="$@"; [ -z "$names" ] && names=$(ls -d */ | tr -d /)
printf "%-28s %6s %10s %10s %6s %6s %s\n" run steps last_conv last_fail negJ retry status
for n in $names; do
  L="$n/$n.log"; [ -f "$L" ] || { printf "%-28s (no log)\n" $n; continue; }
  steps=$(grep -c "^------- converged at time" "$L")
  fails=$(grep -c "failed to converge" "$L")
  lc=$(grep "^------- converged at time" "$L" | tail -1 | awk '{print $NF}')
  lf=$(grep "failed to converge at time" "$L" | tail -1 | awk '{print $NF}')
  nj=$(grep -c "negative jacobian" "$L")
  st="running"
  grep -q "N O R M A L   T E R M" "$L" && st="NORMAL_TERMINATION"
  grep -q "E R R O R   T E R M" "$L" && st="ERROR_TERMINATION"
  printf "%-28s %6s %10s %10s %6s %6s %s\n" $n "$steps" "${lc:--}" "${lf:--}" $nj $fails $st
done
