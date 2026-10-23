"""Batch 416 (2026-10-22; the user on L413b: "I still see a little poking through ... try things like increasing the penalty,
adding an offset distance, etc."). Sets parameters of named contacts, one change per variant (NOT IN SOURCE: these contacts are
ours). Each KEY=VALUE replaces <KEY>...</KEY> inside <contact name="NAME" ...>; a KEY missing there is added before </contact>.
usage: py -3.10 build_batch416.py SRC_RUN NEW_RUN NAME:KEY=VALUE [NAME:KEY=VALUE ...]"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from junction_view_data import load  # noqa: E402
from build_batch162 import write  # noqa: E402


def main(src, new, sets):
    t, X, E = load(src)
    done = []
    for s in sets:
        name, kv = s.split(':', 1); key, val = kv.split('=', 1)
        m = re.search(r'<contact name="%s"[^>]*>.*?</contact>' % re.escape(name), t, re.S)
        assert m, name
        blk = m.group(0)
        old = re.search(r'<%s>([^<]*)</%s>' % (key, key), blk)
        nb = (blk.replace(old.group(0), f'<{key}>{val}</{key}>', 1) if old else
              blk.replace('</contact>', f'\t<{key}>{val}</{key}>\n\t\t</contact>'))
        t = t[:m.start()] + nb + t[m.end():]
        done.append(f'{name}.{key} {old.group(1) if old else "(unset)"} -> {val}')
    msg = 'NOT IN SOURCE (contact settings): ' + '; '.join(done)
    print(new); print(msg)
    write(new, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch416.py)\n{msg}\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
