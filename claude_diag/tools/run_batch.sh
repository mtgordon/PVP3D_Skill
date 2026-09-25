#!/bin/bash
# usage: run_batch.sh THREADS name1 name2 ...   (each name = claude_diag/runs/<name>/<name>.feb)
FEBIO="/c/Program Files/FEBioStudio/bin/febio4.exe"
RUNS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../runs" && pwd)"
T=$1; shift
for n in "$@"; do
  (
    cd "$RUNS/$n" || exit 1
    rm -f "$n.log" "$n.xplt" "$n.stdout"
    OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T "$FEBIO" -i "$n.feb" -silent > "$n.stdout" 2>&1
    echo "$n exit=$?" >> "$RUNS/_finished.txt"
  ) &
done
wait
