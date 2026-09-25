# PVP3D: Abaqus-to-FEBio skill and model-building tools

Material from converting the PVP3D pelvic-floor model (levator ani, vaginal wall, perineal body, and their supports) from
Abaqus to FEBio, and from getting the FEBio version to converge. The FEBio model itself will be added later.

## Contents

| Path | What it is |
|---|---|
| `skill/abaqus-febio-fea-pipeline/` | A Claude Code skill: comparing and converting Abaqus (`.inp`) and FEBio (`.feb`) models, tracing a converted model back to its Abaqus source, and debugging FEBio convergence. `SKILL.md` is the entry point; `references/` holds the detail and `scripts/` the helper tools. |
| `PVP3DModel_job.inp` | The original Abaqus model: `*Dynamic, Explicit` over 1 s, pressures ramped with a smooth step, 291 `CONN3D2` connectors for the ligaments and attachments. |
| `claude_diag/tools/` | Python tools that build the FEBio model variants from a base `.feb` (`variants*.py`, `build_batch*.py`), fit the lofted surfaces' materials to the Abaqus connectors they replace (`loft_survey.py`, `loft_fit_all.py`, `cl_usl_section.py`), run batches (`run_batch.sh`, `autostop_runs.sh`) and compare results (`batch20_table.py`, `tie_gap.py`). |
| `claude_diag/README_2026-09-24b.md` | Notes from the latest round of runs (batches 20-27): what changed, what converged, and what was learned. |

## Using the skill

Copy `skill/abaqus-febio-fea-pipeline/` into your Claude Code skills folder (`~/.claude/skills/`, on Windows
`C:\Users\<you>\.claude\skills\`). Claude Code then offers it for work on `.inp` / `.feb` files. The scripts need
Python 3 with numpy (on Windows, `py -0p` lists the installed interpreters).

## Using the tools

- Python 3.10 with numpy; FEBio 4 (`febio4.exe` from FEBio Studio; `run_batch.sh` expects
  `C:/Program Files/FEBioStudio/bin/febio4.exe`).
- `claude_diag/tools/paths.py` finds everything relative to its own location: the Abaqus source at the repository root,
  and runs in `claude_diag/runs/<name>/<name>.feb`. The build scripts start from base models that will be added with the
  FEBio model.
- Every variant is written with a `<name>.feb.changes.txt` next to it that lists each change, and marks anything not in
  the Abaqus source as NOT IN SOURCE.

## Licence

No licence has been chosen yet, so all rights are reserved by the owner. Ask before reusing or redistributing.
