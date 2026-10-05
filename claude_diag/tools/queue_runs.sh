#!/bin/bash
# usage: queue_runs.sh MAXPROC THREADS name1 name2 ...
# Starts each named run (runs/<name>/<name>.feb, via run_batch.sh THREADS name) once fewer than MAXPROC febio4.exe
# processes are running, in the order given, each with its own 30-min stall monitor (autostop_runs.sh 30 300 name).
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "$TOOLS/../runs" && pwd)"
M=$1; T=$2; shift 2
for n in "$@"; do
  while true; do
    c=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | Measure-Object).Count" | tr -d '\r')
    [ "$c" -lt "$M" ] && break
    sleep 30
  done
  echo "$(date '+%F %T') queue start ($T threads): $n" >> "$RUNS/_autostop.txt"
  nohup bash "$TOOLS/run_batch.sh" "$T" "$n" > /dev/null 2>&1 &
  sleep 3
  nohup bash "$TOOLS/autostop_runs.sh" 30 300 "$n" > /dev/null 2>&1 &
  sleep 20
done
