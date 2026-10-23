#!/bin/bash
# usage: pace_vs.sh BASE_RUN run ...   For each run: converged steps, last converged t, and BASE_RUN's t at the same step count
# (2026-10-14: a run "converging" 1e-9 steps is never stopped by autostop_runs.sh; this shows a crawl against a known run).
RUNS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
B=$1; shift
for r in "$@"; do
  n=$(grep -c "converged at time" "$RUNS/$r/$r.log" 2>/dev/null); t=$(grep "converged at time" "$RUNS/$r/$r.log" 2>/dev/null | tail -1 | awk '{print $NF}')
  tb=$(grep "converged at time" "$RUNS/$B/$B.log" | sed -n "${n:-1}p" | awk '{print $NF}')
  printf '%-32s step %5s  t %-11s | %s at that step: t %s\n' "$r" "$n" "$t" "$B" "$tb"
done
