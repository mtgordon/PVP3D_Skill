"""Batch 249 (2026-10-09, unattended; the range-of-motion goal): the tube line barely moves (L245_tube_cx3_src_*: Ba / Bp /
C motion 10-16 mm vs 30-60 mm on the ledge line). Diagnostics, one change each on a tube run:
  cdamp200   the curves' and the apex cap's mass damping C 1000 -> 200 (= the walls'): is the damping holding them back?
  curvesoft  the curves' wall-material halves (canal_wrap_hex_avw, canal_wrap_hex_pvw; not the cervix band or the cap) on their
             own copy of the wall material at x0.3 (Curve_soft): is the closed tube's stiffness what restricts the motion?
Both NOT IN SOURCE / DIAGNOSTIC.
  curvesoftF (e.g. curvesoft0.6) the same at x F (2026-10-09 night: x0.3 crawled from t 0.67 on negative Jacobians)
usage: py -3.10 build_batch249.py SRC NAME cdamp200|curvesoft|curvesoftF"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants7 import JOBS  # noqa: E402
from build_batch162 import write  # noqa: E402

RUNS = os.path.join(JOBS, 'claude_diag', 'runs')
CURVE_DOMS = ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw', 'canal_wrap_hex_avw_cervix', 'apex_contour')

if __name__ == '__main__':
    src, name, mode = sys.argv[1:4]
    t = open(os.path.join(RUNS, src, src + '.feb'), encoding='utf-8').read()
    log = []
    if mode == 'cdamp200':
        for d in CURVE_DOMS:
            pat = r'(<body_load name="wall_damping_%s"[^>]*>\s*<C>)([^<]+)(</C>)' % d
            m = re.search(pat, t)
            if m:
                t = t[:m.start()] + m.group(1) + '200' + m.group(3) + t[m.end():]
                log.append(f'{d} C {m.group(2)} -> 200')
        msg = 'DIAGNOSTIC: curve / cap mass damping -> 200 (the walls\' value): ' + '; '.join(log)
    elif mode.startswith('curvesoft'):
        F = float(mode[len('curvesoft'):] or 0.3)
        m = re.search(r'\t\t<material id="(\d+)" name="Vagina_AVW" type="Yeoh">.*?</material>\n', t, re.S)
        blk = m.group(0)
        nid = max(int(v) for v in re.findall(r'<material id="(\d+)"', t)) + 1
        new = blk.replace(f'id="{m.group(1)}" name="Vagina_AVW"', f'id="{nid}" name="Curve_soft"')
        for k in ('c1', 'c2', 'k'):
            old = float(re.search(r'<%s>([^<]+)</%s>' % (k, k), new).group(1))
            new = re.sub(r'<%s>[^<]+</%s>' % (k, k), f'<{k}>{old * F:.7g}</{k}>', new)
            log.append(f'{k} {old:.7g} -> {old * F:.7g}')
        i = t.index('\t</Material>')
        t = t[:i] + new + t[i:]
        for d in ('canal_wrap_hex_avw', 'canal_wrap_hex_pvw'):
            t, n = re.subn(r'<SolidDomain name="%s" mat="[^"]+" />' % d, f'<SolidDomain name="{d}" mat="Curve_soft" />', t)
            assert n == 1, d
        msg = ('DIAGNOSTIC, NOT IN SOURCE: canal_wrap_hex_avw / _pvw on Curve_soft = Vagina_AVW x' + f'{F:g}' + ' (' + ', '.join(log)
               + '); the cervix band and the cap unchanged')
    print(msg)
    write(name, t, f'base: claude_diag/runs/{src}/{src}.feb (tools/build_batch249.py)\n{msg}\n')
