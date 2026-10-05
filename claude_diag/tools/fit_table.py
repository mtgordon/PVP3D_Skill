"""The paper-fit sensitivity table (2026-10-02): for every finished run of a batch, Ba / Bp / C (summary_table.py, Bump et al.
1996; + below the hymen) and the change from the reference case (the same case on the reference base), cached in
claude_diag/fit_cache.csv (a run is re-read only if its log is newer than its cache row).
usage: py -3.10 fit_table.py [--batch 150] [--ref-healthy L149_seamspr_pm_outer_fast] [--ref-P1 L150_outerfast_P1]
                             [--ref-P2 L150_outerfast_P2] [RUN ...]"""
import argparse
import csv
import os
import re
import subprocess
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(os.path.dirname(TOOLS), 'runs')
CACHE = os.path.join(os.path.dirname(TOOLS), 'fit_cache.csv')
TARGET = {'healthy': 'Ba, Bp < 0', 'P1': 'Ba < 0, Bp +4', 'P2': 'Ba < 0, Bp +9'}


def load():
    if not os.path.exists(CACHE):
        return {}
    with open(CACHE, encoding='utf-8', newline='') as f:
        return {r['run']: r for r in csv.DictReader(f)}


def save(c):
    with open(CACHE, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['run', 'mtime', 't', 'failed', 'wall', 'Ba', 'Bp', 'C'])
        w.writeheader()
        w.writerows(c.values())


def measure(runs, cache):
    todo = []
    for r in runs:
        log = os.path.join(RUNS, r, r + '.log')
        if not os.path.exists(log) or 'T E R M' not in open(log, errors='ignore').read()[-3000:]:
            continue
        mt = str(int(os.path.getmtime(log)))
        if r not in cache or cache[r]['mtime'] != mt:
            todo.append((r, mt))
    if todo:
        out = subprocess.run([sys.executable, os.path.join(TOOLS, 'summary_table.py')] + [r for r, _ in todo],
                             capture_output=True, text=True).stdout
        for r, mt in todo:
            m = re.search(r'^\| `%s` \| ([^|]+)\| ([^|]+)\|[^|]+\| ([^|]+)\|.*\| ([^|]+)\| ([^|]+)\| ([^|]+)\|\s*$' % re.escape(r),
                          out, re.M)
            if m:
                t, fl, wall, ba, bp, c = (g.strip() for g in m.groups())
                cache[r] = dict(run=r, mtime=mt, t=t, failed=fl, wall=wall, Ba=ba, Bp=bp, C=c)
        save(cache)
    return cache


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--batch', default='150')
    ap.add_argument('--ref-healthy', default='L149_seamspr_pm_outer_fast')
    ap.add_argument('--ref-P1', default='L150_outerfast_P1')
    ap.add_argument('--ref-P2', default='L150_outerfast_P2')
    ap.add_argument('runs', nargs='*')
    a = ap.parse_args()
    ref = {'healthy': a.ref_healthy, 'P1': a.ref_P1, 'P2': a.ref_P2}
    runs = a.runs or sorted(d for d in os.listdir(RUNS) if re.match(r'L%s_\w+x(05|2)_(healthy|P1|P2)$' % a.batch, d))
    cache = measure(list(ref.values()) + runs, load())
    print('| run | case | t | failed | wall | Ba / Bp / C | dBa / dBp / dC from the reference | target |')
    print('|---|---|---|---|---|---|---|---|')
    for r in list(ref.values()) + runs:
        if r not in cache:
            continue
        case = next(c for c in ('healthy', 'P1', 'P2') if r.endswith(c) or r == ref[c])
        row = cache[r]
        d = ''
        if r != ref[case] and ref[case] in cache:
            d = ' / '.join(f'{float(row[k]) - float(cache[ref[case]][k]):+.1f}' for k in ('Ba', 'Bp', 'C'))
        print(f'| `{r}` | {case} | {row["t"]} | {row["failed"]} | {row["wall"]} | {row["Ba"]} / {row["Bp"]} / {row["C"]} | '
              f'{d} | {TARGET[case]} |')
