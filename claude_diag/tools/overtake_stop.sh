#!/bin/bash
# usage: overtake_stop.sh NEW OLD [MAXCHECKS] [INTERVAL_S]  -- every INTERVAL_S (default 600) s compare the last converged time of runs/NEW and runs/OLD; when NEW has
# passed OLD while OLD is still running, stop OLD (stop_runs.sh) and note it in runs/_autostop.txt. Ends when either finishes.
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
NEW=$1; OLD=$2; MAX=${3:-11}; DT=${4:-600}
lastt() { grep -o 'converged at time *: *[0-9.e+-]*' "$RUNS/$1/$1.log" 2>/dev/null | tail -1 | awk -F: '{print $2+0}'; }
done_() { tail -c 400 "$RUNS/$1/$1.log" 2>/dev/null | grep -q "T E R M I N A T I O N"; }
for ((i = 0; i < MAX; i++)); do
  if done_ "$OLD"; then echo "$OLD finished first"; exit 0; fi
  if done_ "$NEW"; then echo "$NEW finished"; exit 0; fi
  tn=$(lastt "$NEW"); to=$(lastt "$OLD")
  if awk -v a="$tn" -v b="$to" 'BEGIN{exit !(a > b)}'; then
    bash "$TOOLS/stop_runs.sh" "$OLD"
    echo "$(date '+%F %T') stopped $OLD at t $to: $NEW passed it (t $tn) (the user: kill the old one once the new one passes it)" >> "$RUNS/_autostop.txt"
    echo "stopped $OLD at t $to (new at $tn)"; exit 0
  fi
  sleep "$DT"
done
echo "no overtake within $MAX checks (new $(lastt "$NEW"), old $(lastt "$OLD"))"
