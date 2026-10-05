"""Add or update rows of claude_diag/run_log_notes.csv (the hand-written part of PVP3D_run_log.xlsx; run_log.py rebuilds it).
usage: py -3.10 run_notes.py RUN [--why TEXT] [--result TEXT] [--label TEXT] [--focus rotation|pvw|la|] [--source TEXT]
A new run gets a row; for an existing run only the given fields change. The file keeps its UTF-8 BOM and column order."""
import csv
import os
import sys

NOTES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'run_log_notes.csv')
FIELDS = ['run', 'why', 'result', 'label', 'focus', 'source']


def main(argv):
    run, rest = argv[0], argv[1:]
    upd = {}
    while rest:
        k, v, rest = rest[0], rest[1], rest[2:]
        assert k.startswith('--') and k[2:] in FIELDS[1:], k
        upd[k[2:]] = v
    with open(NOTES, encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    row = next((r for r in rows if r['run'] == run), None)
    if row is None:
        row = {k: '' for k in FIELDS}
        row['run'] = run
        rows.append(row)
        print('added', run)
    else:
        print('updated', run, sorted(upd))
    row.update(upd)
    with open(NOTES, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main(sys.argv[1:])
