#!/bin/bash
# usage: wait_slot.sh MAXPROC MINUTES [name ...]
# A capped wait (the user's rule, 2026-09-27: never wait more than ~30 min without re-checking every run). Returns as
# soon as fewer than MAXPROC febio4.exe processes run (a run ended or was stopped: a slot is free), or after MINUTES,
# then prints progress.sh for the given runs (default: the runs whose febio4.exe was running at the start).
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M=$1; MIN=$2; shift 2
running() {
  powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | ForEach-Object { \$_.CommandLine }" 2>/dev/null \
    | tr -d '\r' | grep -o '[A-Za-z0-9_]*\.feb' | sed 's/\.feb$//'
}
names="$*"; [ -z "$names" ] && names=$(running | tr '\n' ' ')
end=$((SECONDS + MIN * 60)); why="cap of $MIN min reached"
while [ $SECONDS -lt $end ]; do
  c=$(running | wc -l)
  if [ "$c" -lt "$M" ]; then why="only $c febio4.exe running (< $M)"; break; fi
  sleep 60
done
echo "$(date '+%F %T') wait over: $why"
bash "$TOOLS/progress.sh" $names
