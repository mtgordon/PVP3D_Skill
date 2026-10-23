#!/bin/bash
# usage: wait_t.sh MINUTES T run ...   Returns when any run's last converged t >= T, or any run ends / is gone (wait_any_end.sh
# rules), or after MINUTES (the 30-min cap), then prints each run's last converged t.
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; RUNS="$(cd "$TOOLS/../runs" && pwd)"
M=$1; T=$2; shift 2; end=$(( $(date +%s) + M*60 ))
lt(){ grep "converged at time : [0-9]" "$RUNS/$1/$1.log" 2>/dev/null | tail -1 | awk '{print $NF}'; }
while [ $(date +%s) -lt $end ]; do
  for r in "$@"; do
    t=$(lt $r); awk -v a="${t:-0}" -v b="$T" 'BEGIN{exit !(a>=b)}' && { echo "$r reached t $t"; break 2; }
    grep -q "T E R M" "$RUNS/$r/$r.log" 2>/dev/null && { echo "$r ended"; break 2; }
  done
  sleep 60
done
for r in "$@"; do echo "$r $(lt $r) $(grep -o 'N O R M A L\|E R R O R' "$RUNS/$r/$r.log" 2>/dev/null | head -1)"; done
