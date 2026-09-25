# PVP3D: Abaqus-to-FEBio skill and model-building tools

Material from converting the PVP3D pelvic-floor model (levator ani, vaginal wall, perineal body, and their supports) from
Abaqus to FEBio, and from getting the FEBio version to converge, including the FEBio model that runs most reliably.

## Contents

| Path | What it is |
|---|---|
| `skill/abaqus-febio-fea-pipeline/` | A Claude Code skill: comparing and converting Abaqus (`.inp`) and FEBio (`.feb`) models, tracing a converted model back to its Abaqus source, and debugging FEBio convergence. `SKILL.md` is the entry point; `references/` holds the detail and `scripts/` the helper tools. |
| `PVP3DModel_job.inp` | The original Abaqus model: `*Dynamic, Explicit` over 1 s, pressures ramped with a smooth step, 291 `CONN3D2` connectors for the ligaments and attachments. |
| `claude_diag/tools/` | Python tools that build the FEBio model variants from a base `.feb` (`variants*.py`, `build_batch*.py`), fit the lofted surfaces' materials to the Abaqus connectors they replace (`loft_survey.py`, `loft_fit_all.py`, `cl_usl_section.py`), run batches (`run_batch.sh`, `autostop_runs.sh`) and compare results (`batch20_table.py`, `tie_gap.py`). |
| `claude_diag/runs/L12_stabdamp_p3/` | **The FEBio model** (`L12_stabdamp_p3.feb`), its log, its change list, and `L12_stabdamp_p3_history.md`: every change from the converted file to this one. See below. |
| `claude_diag/README_2026-09-24b.md` | Notes from the latest round of runs (batches 20-27): what changed, what converged, and what was learned. |

## The FEBio model: `L12_stabdamp_p3`

The most stable version so far: it runs to full load (t = 1) with **no failed time steps**, in about 12 minutes on
10 threads (`febio4 -i L12_stabdamp_p3.feb`; FEBio 4, implicit dynamic analysis). It is the Abaqus model converted to
FEBio, with:

- the connector families replaced by lofted shell surfaces ("lofts"), each with its own 1-term Ogden material fitted
  to the Abaqus connectors it replaces (strip model); the lofts standing in for near-zero connectors are very soft and
  at tissue density;
- the other lofts still carry the density inherited from the Abaqus beam section they were cloned from (7.8e-07
  t/mm^3, ~1.35 kg in all, while the connectors are massless). Tissue density gave the same solution but converged
  far worse, so it was left for now;
- the whole anchored edge of the PM_PeB lofts fixed (a loft convention);
- the line search on (`lstol` 0.9), which gave 16x fewer failed steps than without it;
- the stabilizing ground springs of the early conversion removed (the Abaqus model has none);
- **NOT IN SOURCE:** every surface pressure at 1/3 of the Abaqus value (0.00467 instead of 0.014 MPa), while the
  model's deformation is being checked;
- **NOT IN SOURCE:** mass damping (C = 20/s) from the start, to keep the freed vaginal wall from moving too fast for the
  solver;
- the posterior-arcus spring chain with 10x the Abaqus truss mass (a screening shortcut; the Abaqus chain mass is
  the like-for-like version).

At t = 1 the levator ani displacement is 6.3 / 10.2 / 13.5 mm (median / p90 / max). A 1 s dynamic ramp can lag the
load (see the notes), so whether this is the equilibrium is being checked with the load held to t = 2. The plot file
(`.xplt`, ~170 MB) is over GitHub's size limit: run the model to regenerate it.

## Using the skill

Copy `skill/abaqus-febio-fea-pipeline/` into your Claude Code skills folder (`~/.claude/skills/`, on Windows
`C:\Users\<you>\.claude\skills\`). Claude Code then offers it for work on `.inp` / `.feb` files. The scripts need
Python 3 with numpy (on Windows, `py -0p` lists the installed interpreters).

## Using the tools

- Python 3.10 with numpy; FEBio 4 (`febio4.exe` from FEBio Studio; `run_batch.sh` expects
  `C:/Program Files/FEBioStudio/bin/febio4.exe`).
- `claude_diag/tools/paths.py` finds everything relative to its own location: the Abaqus source at the repository root,
  and runs in `claude_diag/runs/<name>/<name>.feb`. The build scripts start from the earlier models in the chain (see
  the history file), which are not included; `L12_stabdamp_p3.feb` can serve as the base for new variants.
- Every variant is written with a `<name>.feb.changes.txt` next to it that lists each change, and marks anything not in
  the Abaqus source as NOT IN SOURCE.

## Licence

No licence has been chosen yet, so all rights are reserved by the owner. Ask before reusing or redistributing.
