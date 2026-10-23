"""One-off (2026-09-24): point the F: copy of the next-session prompt at the new project location."""
import os
import re

OLD = 'C:\\Users\\mgordon\\Downloads\\jobs'
NEW = 'F:\\Developing PVP3D Model for FEBio'
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
p = os.path.join(root, 'NEXT_SESSION_PROMPT_2026-09-23c.md')
s = open(p, 'rb').read().decode('utf-8')
n = s.count(OLD)
s = s.replace(OLD, NEW)
note = ('\n**The project moved on 2026-09-24** from `C:\\Users\\mgordon\\Downloads\\jobs` to `F:\\Developing PVP3D Model for FEBio`\n'
        '(an external drive). The Abaqus source is copied to `PVP3DModel_job.inp` in that folder. The tools now find their\n'
        'paths from their own location (`claude_diag/tools/paths.py`), so they work wherever the folder lives. Use this F: folder.\n'
        'The Box results folder is unchanged.\n')
s = s.replace('**Step 1: invoke the correct skill', note.lstrip('\n') + '\n**Step 1: invoke the correct skill', 1)
open(p, 'wb').write(s.encode('utf-8'))
print(f'{p}: {n} path(s) updated, move note added')
