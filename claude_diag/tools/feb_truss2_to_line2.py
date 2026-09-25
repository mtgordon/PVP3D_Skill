"""Make a .feb readable by FEBio Studio 3.2: rewrite truss elements from the `truss2` + <SolidDomain> form (accepted
by the FEBio 4.13 solver, rejected by FEBio Studio: 'tag "Elements" : invalid value for attribute "type"') to the
documented FEBio 4 form, `line2` + <BeamDomain type="linear-truss">. Only those lines change; every other byte of
the file is kept. The solver gives identical results for both forms (claude_diag/mini/truss_mass.feb vs
truss_mass_line2.feb).

usage: py -3.10 feb_truss2_to_line2.py IN.feb [OUT.feb]      (default OUT: IN_FEBioStudio.feb; never overwrites)
"""
import os
import re
import sys

src = sys.argv[1]
dst = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(src)[0] + '_FEBioStudio.feb'
if os.path.exists(dst):
    sys.exit(f'{dst} exists; not overwriting')
text = open(src, 'rb').read().decode('latin-1')
names = re.findall(r'<Elements type="truss2" name="([^"]+)">', text)
if not names:
    sys.exit('no truss2 element blocks; nothing to do')
n_el = len(names)
text = re.sub(r'<Elements type="truss2" name="', '<Elements type="line2" name="', text)
n_dom = 0
for nm in names:
    pat = re.compile(r'<SolidDomain name="' + re.escape(nm) + r'" mat="([^"]+)">(.*?)</SolidDomain>', re.S)
    m = pat.search(text)
    assert m and 'cross_sectional_area' in m.group(2), f'domain for {nm} not found or has no cross_sectional_area'
    text = text[:m.start()] + f'<BeamDomain name="{nm}" mat="{m.group(1)}" type="linear-truss">{m.group(2)}</BeamDomain>' \
        + text[m.end():]
    n_dom += 1
open(dst, 'wb').write(text.encode('latin-1'))
print(f'wrote {dst}: {n_el} Elements block(s) truss2 -> line2, {n_dom} SolidDomain -> BeamDomain type="linear-truss" ({names})')
