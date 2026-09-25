#!/bin/bash
# usage: watch_runs.sh name1 name2 ...  -- prints one line per milestone / termination / stall
RUNS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
declare -A mile stall lastt lastchg done
marks="0.1 0.2 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0"
now(){ date +%s; }
first=1
for n in "$@"; do mile[$n]=0; stall[$n]=0; lastt[$n]=""; lastchg[$n]=$(now); done
while true; do
  alive=0
  for n in "$@"; do
    [ -n "${done[$n]}" ] && continue
    alive=1
    L="$RUNS/$n/$n.log"; [ -f "$L" ] || continue
    t=$(grep "^------- converged at time" "$L" | tail -1 | awk '{print $NF}')
    r=$(grep -c "failed to converge" "$L")
    if [ "$t" != "${lastt[$n]}" ]; then lastt[$n]=$t; lastchg[$n]=$(now); stall[$n]=0; fi
    for m in $marks; do
      if awk -v a="${t:-0}" -v b=$m 'BEGIN{exit !(a>=b)}' && awk -v a="${mile[$n]}" -v b=$m 'BEGIN{exit !(a<b)}'; then
        mile[$n]=$m; [ -n "$first" ] || echo "$(date +%H:%M) $n passed t=$m (t=$t, retries=$r)"
      fi
    done
    if grep -q "N O R M A L   T E R M" "$L"; then echo "$(date +%H:%M) $n NORMAL TERMINATION t=$t retries=$r"; done[$n]=1; fi
    if grep -q "E R R O R   T E R M" "$L"; then echo "$(date +%H:%M) $n ERROR TERMINATION last converged t=$t retries=$r"; done[$n]=1; fi
    if [ -z "${done[$n]}" ] && [ ${stall[$n]} -eq 0 ] && [ $(( $(now) - ${lastchg[$n]} )) -gt 1200 ]; then
      echo "$(date +%H:%M) $n STALLED: no new converged step for 20 min at t=$t retries=$r"; stall[$n]=1
    fi
  done
  [ $alive -eq 0 ] && exit 0
  first=""
  sleep 60
done
