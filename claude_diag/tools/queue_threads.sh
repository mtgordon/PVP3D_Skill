#!/bin/bash
# usage: queue_threads.sh LIMIT THREADS name1 name2 ...
# Like queue_runs.sh, but counts THREADS, not processes (2026-10-09: 12-thread runs still going while an 8-thread batch is
# queued). Each run's thread count is read from runs/_threads.txt ("name threads", appended here at start; runs not listed
# there count as 12, queue_runs.sh's usual). A run starts once the running total + THREADS <= LIMIT (<= 64, the user's rule).
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
LIM=$1; T=$2; shift 2
[ "$LIM" -gt 64 ] && { echo "refused: LIMIT $LIM > 64" >&2; exit 1; }
LOCK="$RUNS/.queue_lock"
total() {
  powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | ForEach-Object { \$_.CommandLine }" 2>/dev/null \
    | tr -d '\r' | grep -o '[A-Za-z0-9_]*\.feb' | sed 's/\.feb$//' | while read n; do
      t=$(grep "^$n " "$RUNS/_threads.txt" 2>/dev/null | tail -1 | cut -d' ' -f2); echo "${t:-12}"; done | awk '{s+=$1} END {print s+0}'
}
for n in "$@"; do
  while true; do
    if mkdir "$LOCK" 2>/dev/null; then
      [ $(( $(total) + T )) -le "$LIM" ] && break
      rmdir "$LOCK"
    elif [ -n "$(find "$LOCK" -maxdepth 0 -mmin +10 2>/dev/null)" ]; then
      rmdir "$LOCK" 2>/dev/null
    fi
    sleep 30
  done
  echo "$n $T" >> "$RUNS/_threads.txt"
  echo "$(date '+%F %T') queue start ($T threads): $n" >> "$RUNS/_autostop.txt"
  nohup bash "$TOOLS/run_batch.sh" "$T" "$n" > /dev/null 2>&1 &
  sleep 3
  nohup bash "$TOOLS/autostop_runs.sh" 30 300 "$n" > /dev/null 2>&1 &
  sleep 20
  rmdir "$LOCK"
done
