#!/bin/bash
# usage: autostop_runs.sh STALL_MIN MAX_MIN name1 name2 ...
# Like wait_runs.sh, but applies the user's rule (2026-09-24): a run with no new converged step for STALL_MIN minutes
# is STOPPED (stop_runs.sh; its files are kept), so it frees its cores for the others. Exits once every run has
# terminated, exited, or been stopped (or after MAX_MIN, stopping nothing more), then prints one line per run.
# Events go to runs/_autostop.txt.
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
S=$1; M=$2; shift 2
declare -A lastt lastchg state
now(){ date +%s; }
t0=$(now)
for n in "$@"; do lastt[$n]=""; lastchg[$n]=$(now); state[$n]=""; done
echo "$(date '+%F %T') autostop start: stall ${S} min, max ${M} min: $*" >> "$RUNS/_autostop.txt"
while true; do
  procs=$(powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | ForEach-Object { \$_.CommandLine }" 2>/dev/null)
  pending=0
  for n in "$@"; do
    [ -n "${state[$n]}" ] && continue
    L="$RUNS/$n/$n.log"
    t=$(grep "^------- converged at time" "$L" 2>/dev/null | tail -1 | awk '{print $NF}')
    if [ "$t" != "${lastt[$n]}" ]; then lastt[$n]=$t; lastchg[$n]=$(now); fi
    if grep -q "N O R M A L   T E R M" "$L" 2>/dev/null; then state[$n]="finished"; continue; fi
    if grep -q "E R R O R   T E R M" "$L" 2>/dev/null; then state[$n]="failed (error termination)"; continue; fi
    if ! grep -q "$n\.feb" <<< "$procs"; then state[$n]="exited without a termination line"; continue; fi
    if [ $(( $(now) - ${lastchg[$n]} )) -gt $((S * 60)) ]; then
      bash "$TOOLS/stop_runs.sh" "$n" >> "$RUNS/_autostop.txt" 2>&1
      state[$n]="stalled ${S}+ min, stopped $(date '+%H:%M')"
      echo "$(date '+%F %T') $n stopped: no new converged step since t=${t:-0}" >> "$RUNS/_autostop.txt"
      continue
    fi
    pending=1
  done
  [ $pending -eq 0 ] && break
  [ $(( $(now) - t0 )) -gt $((M * 60)) ] && { echo "MAX_MIN reached"; break; }
  sleep 60
done
for n in "$@"; do
  L="$RUNS/$n/$n.log"
  t=$(grep "^------- converged at time" "$L" 2>/dev/null | tail -1 | awk '{print $NF}')
  r=$(grep -c "failed to converge" "$L" 2>/dev/null)
  echo "$n t=${t:-0} retries=${r:-0} ${state[$n]:-running}"
done | tee -a "$RUNS/_autostop.txt"
