#!/bin/bash
# usage: queue_runs.sh MAXPROC THREADS name1 name2 ...
# Starts each named run (runs/<name>/<name>.feb, via run_batch.sh THREADS name) once fewer than MAXPROC febio4.exe
# processes are running, in the order given, each with its own 30-min stall monitor (autostop_runs.sh 30 300 name).
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
M=$1; T=$2; shift 2
# the user, 2026-10-06: never oversubscribe the machine (64 logical cores): MAXPROC x THREADS must stay <= 64, counting
# every febio4.exe already running (use the same THREADS for everything queued, or lower MAXPROC to fit)
if [ $((M * T)) -gt 64 ]; then echo "refused: MAXPROC $M x THREADS $T = $((M * T)) > 64 logical cores" >&2; exit 1; fi
# 2026-10-07: two queues started runs in the same second (6 x 12 threads). The count-and-start step now holds a lock
# (runs/.queue_lock, a directory; a lock older than 10 min is stale and removed) shared by every queue.
LOCK="$RUNS/.queue_lock"
for n in "$@"; do
  while true; do
    if mkdir "$LOCK" 2>/dev/null; then
      c=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | Measure-Object).Count" | tr -d '\r')
      [ "$c" -lt "$M" ] && break
      rmdir "$LOCK"
    elif [ -n "$(find "$LOCK" -maxdepth 0 -mmin +10 2>/dev/null)" ]; then
      rmdir "$LOCK" 2>/dev/null
    fi
    sleep 30
  done
  echo "$(date '+%F %T') queue start ($T threads): $n" >> "$RUNS/_autostop.txt"
  nohup bash "$TOOLS/run_batch.sh" "$T" "$n" > /dev/null 2>&1 &
  sleep 3
  nohup bash "$TOOLS/autostop_runs.sh" 30 300 "$n" > /dev/null 2>&1 &
  sleep 20          # febio4.exe is up by now, so the next count sees it
  rmdir "$LOCK"
done
