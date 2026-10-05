#!/bin/bash
# usage: wait_any_end.sh MINUTES name ...   Returns when any named run's log gets a termination line or its febio4.exe
# is gone (ended, stopped or crashed), or after MINUTES (the 30-min cap rule), then prints progress.sh for them.
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
MIN=$1; shift
state() {
  procs=$(powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | ForEach-Object { \$_.CommandLine }" 2>/dev/null)
  for n in "$@"; do
    L="$RUNS/$n/$n.log"
    if grep -q "T E R M" "$L" 2>/dev/null; then echo "$n:end"; elif [ -f "$L" ] && ! grep -q "$n\.feb" <<< "$procs"; then echo "$n:gone"; fi
  done | sort | tr '\n' ' '
}
s0=$(state "$@"); end=$((SECONDS + MIN * 60)); why="cap of $MIN min reached"
while [ $SECONDS -lt $end ]; do
  sleep 60
  s=$(state "$@")
  if [ "$s" != "$s0" ]; then why="changed: $s"; break; fi
done
echo "$(date '+%F %T') wait over: $why"
bash "$TOOLS/progress.sh" "$@"
