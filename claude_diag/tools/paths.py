"""Project paths, derived from this file's location, so the project folder can live anywhere (it moved from
one folder to another on 2026-09-24).
tools/ -> claude_diag/ -> project root."""
import os

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
DIAG_DIR = os.path.dirname(TOOLS_DIR)
JOBS_DIR = os.path.dirname(DIAG_DIR)
RUNS_DIR = os.path.join(DIAG_DIR, 'runs')
INP_FILE = os.path.join(JOBS_DIR, 'PVP3DModel_job.inp')     # the Abaqus source
