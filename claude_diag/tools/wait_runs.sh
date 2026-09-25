#!/bin/bash
# usage: wait_runs.sh STALL_MIN MAX_MIN name1 name2 ...
# Exits once every run has terminated or made no new converged step for STALL_MIN minutes (or after MAX_MIN),
# then prints one line per run. One notification instead of a stream of milestones.
RUNS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
S=$1; M=$2; shift 2
declare -A lastt lastchg
now(){ date +%s; }
t0=$(now)
for n in "$@"; do lastt[$n]=""; lastchg[$n]=$(now); done
while true; do
  pending=0
  for n in "$@"; do
    L="$RUNS/$n/$n.log"
    t=$(grep "^------- converged at time" "$L" 2>/dev/null | tail -1 | awk '{print $NF}')
    if [ "$t" != "${lastt[$n]}" ]; then lastt[$n]=$t; lastchg[$n]=$(now); fi
    grep -qE "N O R M A L   T E R M|E R R O R   T E R M" "$L" 2>/dev/null && continue
    [ $(( $(now) - ${lastchg[$n]} )) -gt $((S * 60)) ] && continue
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
  s="running"; [ $(( $(now) - ${lastchg[$n]} )) -gt $((S * 60)) ] && s="stalled ${S}+ min"
  grep -q "N O R M A L   T E R M" "$L" 2>/dev/null && s="finished"
  grep -q "E R R O R   T E R M" "$L" 2>/dev/null && s="failed"
  echo "$n t=${t:-0} retries=${r:-0} $s"
done
