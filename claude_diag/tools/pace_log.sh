#!/bin/bash
# usage: pace_log.sh [EVERY_MIN] [HOURS]   (default 5 min, 10 h; run with nohup ... &)
# Every EVERY_MIN minutes, appends one line per running febio4.exe run to runs/_pace.csv:
#   clock,run,t (last converged),failed attempts,steps
# so crawling runs can be compared by the wall time they take between matching t values (the logs have no clock).
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
E=${1:-5}; HRS=${2:-10}
OUT="$RUNS/_pace.csv"
[ -f "$OUT" ] || echo "clock,run,t,failed,steps" > "$OUT"
end=$((SECONDS + HRS * 3600))
while [ $SECONDS -lt $end ]; do
  now=$(date '+%F %T')
  for n in $(powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | ForEach-Object { \$_.CommandLine }" 2>/dev/null \
             | tr -d '\r' | grep -o '[A-Za-z0-9_]*\.feb' | sed 's/\.feb$//'); do
    L="$RUNS/$n/$n.log"
    t=$(grep "^------- converged at time" "$L" 2>/dev/null | tail -1 | awk '{print $NF}')
    f=$(grep -c "failed to converge" "$L" 2>/dev/null)
    s=$(grep -c "^------- converged at time" "$L" 2>/dev/null)
    echo "$now,$n,${t:-0},$f,$s" >> "$OUT"
  done
  sleep $((E * 60))
done
