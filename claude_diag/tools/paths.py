"""Project paths, derived from this file's location, so the project folder can live anywhere (it moved from
C:/Users/mgordon/Downloads/jobs to F:/Developing PVP3D Model for FEBio on 2026-09-24).
tools/ -> claude_diag/ -> project root."""
import os

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
DIAG_DIR = os.path.dirname(TOOLS_DIR)
JOBS_DIR = os.path.dirname(DIAG_DIR)
RUNS_DIR = os.path.join(DIAG_DIR, 'runs')
INP_FILE = os.path.join(JOBS_DIR, 'PVP3DModel_job.inp')     # copy of the Abaqus source (was in Downloads)

# The skill's helper scripts (xplt_reader, feb_model, ...): $PVP3D_SKILL_SCRIPTS, else the shared repo's .claude/skills/ copy,
# else the user's installed skill (2026-10-23, for sharing the project).
SKILL_SCRIPTS = os.environ.get('PVP3D_SKILL_SCRIPTS') or next(
    (p for p in (os.path.join(JOBS_DIR, '.claude', 'skills', 'abaqus-febio-fea-pipeline', 'scripts'),
                 os.path.join(os.path.expanduser('~'), '.claude', 'skills', 'abaqus-febio-fea-pipeline', 'scripts'))
     if os.path.isdir(p)), os.path.join(os.path.expanduser('~'), '.claude', 'skills', 'abaqus-febio-fea-pipeline', 'scripts'))
FEBIO = os.environ.get('FEBIO4', r'C:\Program Files\FEBioStudio\bin\febio4.exe')
