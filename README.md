# PVP3D: Abaqus-to-FEBio skill and model-building tools

Material from converting the PVP3D pelvic-floor model (levator ani, vaginal wall, perineal body, and their supports) from
Abaqus to FEBio, and from getting the FEBio version to converge, including the current FEBio models, their range of motion
against the paper's prolapse cases, and the rebuilt tube junction. Updated 2026-10-12.

## Contents

| Path | What it is |
|---|---|
| `skill/abaqus-febio-fea-pipeline/` | A Claude Code skill: comparing and converting Abaqus (`.inp`) and FEBio (`.feb`) models, tracing a converted model back to its Abaqus source, and debugging FEBio convergence. `SKILL.md` is the entry point; `references/` holds the detail and `scripts/` the helper tools. |
| `PVP3DModel_job.inp` | The original Abaqus model: `*Dynamic, Explicit` over 1 s, pressures ramped with a smooth step, 291 `CONN3D2` connectors for the ligaments and attachments. |
| `claude_diag/tools/` | Python tools that build the FEBio model variants from a base `.feb` (`variants*.py`, `build_batch*.py`), run batches (`queue_threads.sh`, `wait_any_end.sh`, `autostop_runs.sh`, `stop_runs.sh`) and read results (`rom_report.py`: range of motion; `pop_q_bump.py`: POP-Q points as Bump et al. 1996 define them; `curve_swell.py`: element swelling / folding; `junction_slices.py`: depth-slice pictures of the tube junction). |
| `claude_diag/runs/L258_tubecand_*/` | **The tube line's main model** (2026-10-10): the vaginal canal as a closed tube (walls joined by curved side bands). See below. |
| `claude_diag/runs/L283_junction_fine2z18_*/` | **The tube with the rebuilt junction** (2026-10-12; candidate, not yet adopted). See below. |
| `claude_diag/runs/L228_ext5sl01_pmc_pmd_sph025/`, `claude_diag/runs/L271_ledge_cx3_chains_sheet30_*/` | **The ledge line**: its base, and the version without compressive springs that meets the range-of-motion limits. See below. |
| `claude_diag/README_2026-10-09.md` | Notes for 2026-10-09 to 10-12 (range-of-motion goal, both lines, the junction rebuild: sections 3-7). `README_2026-10-06.md` and `README_2026-10-02.md` before it. |
| `claude_diag/junction_view/` | Depth slices of the tube junction before and after the rebuild. |
| `claude_diag/runs/L149_*`, `L159_*`, `L12_stabdamp_p3/` | Earlier bases and the paper fit (2026-09-24 to 10-05); superseded, kept for history (sections further down). |

## The goal since 2026-10-09: range of motion

The target is Luo et al. 2015 (J Biomech): impairing the supports must be able to move the walls at least as far as the
paper's cases did (Ba >= +10 mm, Bp >= +10 mm below the hymen, apex motion >= 35 mm), while the healthy model stays
unprolapsed (Ba, Bp < 0). Cases: healthy, the paper's P1 / P2 / F6-F8 cases, and MAX (levator 60 %, apical 80 %,
anterior 50 %, posterior 95 % impaired). Shared settings: the remeshed perineal membrane, the cervix x3 (NOT FAITHFUL),
the sphincter-side connectors x0.25 (NOT IN SOURCE), no compression in any connector spring.

## The tube line

`L258_tubecand_*` (main): the walls joined into a closed tube by curved side bands ("curves", x0.6 stiffness), an apex cap,
no pressure on the curves. Healthy Ba / Bp -14.0 / -17.5; MAX +18.4 / +17.2, apex motion 36 mm (every limit met). But in
MAX the curves' ends next to the perineal body swell to 6x their volume and two elements fold inside out: the ends were
point-tied inside the wall's hexes where the wall had been flared into them.

`L283_junction_fine2z18_*` (candidate): the same model with the junction rebuilt (NOT IN SOURCE: `build_batch278.py full`,
then `build_batch280.py --zmax 18 --split 2`): the wall's flare removed (its last layer straight, tied onto the perineal
body's top face), the curves joined to the straight wall and their ends tied onto the perineal body's top face (no point
ties), and every curve ring within 18 mm of the perineal body split in two (the wall rolls over its attachment there in
MAX, which crushed the coarse rings).

| run | Ba / Bp, apex motion | curve swelling, folded elements |
|---|---|---|
| `L258_tubecand_MAX` | +18.4 / +17.2, 36.3 mm | 6.1x, 2 |
| `L283_junction_fine2z18_MAX` | +17.5 / +14.6, 35.5 mm | 3.3x, none |
| `L283_junction_fine2z18_healthy` | -14.7 / -19.1 (main -14.0 / -17.5) | 1.8x, none |

L283 runs ~1.5-2x slower than L258 (its logs, 25-30 MB, are not included). Pictures: `claude_diag/junction_view/`.

## The ledge line

`L228_ext5sl01_pmc_pmd_sph025` is the base (the walls' side seams with short extensions, "ledges", as a convergence aid).
`L271_ledge_cx3_chains_sheet30_*` adds the cervix x3 and turns every connector family into chains of springs inside a
thin sheet at 30 % of its fibre stiffness (NOT IN SOURCE), which holds the healthy wall without compressive springs:
healthy -9.7 / -20.4; MAX +30.7 / +23.9, apex motion 54 mm, smooth to full load. The single-compartment cases fall a
little short (README_2026-10-09.md section 6; 25 % sheets reach a bit further with less healthy margin).

## Earlier: the base of 2026-10-02, `L149_seamspr_pm_outer`

The Abaqus model in FEBio 4 with every connector as a spring (no lofts), run to full load (t = 1). Compared with the
source (everything not in the Abaqus model is marked NOT IN SOURCE in `L149_seamspr_pm_outer.feb.changes.txt`):

- the vaginal walls and perineal body refitted from the source's Marlow data as Yeoh (I1-only, like Marlow);
- the levator at the source's healthy law (PCM-LA_Yamada100%), the levator load at 1/3 of the source;
- **NOT IN SOURCE:** weak springs (0.01 N/mm) across the canal's side seam, where the anterior and posterior walls meet
  edge to edge with only a frictionless contact between them. The source runs this freely as an explicit dynamic
  analysis; FEBio's implicit solver stalls there without them;
- **NOT IN SOURCE:** the perineal membrane (PM_Plane, a display-only part in the source) as a deformable shell clamped
  on its outer arc, with the PM-family connectors tied to it and the 27 PM_conn-family springs anchored on the clamped
  outer arc;
- **NOT IN SOURCE:** wall contacts against the levator (`walls_LA`) and the membrane (`walls_PM`).

## Earlier: the fit to the paper's prolapse cases, `L159_C3full_parcus042_*`

The target is Luo et al. 2015 (J Biomech, rectocele and cystocele): cases P1 and P2, given as POP-Q points. The fit
changes only the healthy model's connective tissue and wall stiffness, within x0.2 to x5. The best set so far, C3,
is the posterior support (P-arcus springs) at x0.42, **NOT FAITHFUL** to the source:

| run (`claude_diag/runs/<run>/<run>.feb`) | Ba / Bp / C (mm; + below the hymen) | paper |
|---|---|---|
| `L159_C3full_parcus042_healthy` | -6.1 / -12.9 / -37.3 | no prolapse |
| `L159_C3full_parcus042_P1` | +9.3 / +5.0 / -25.2 | Ba < 0, Bp +4 |
| `L159_C3full_parcus042_P2` | +14.8 / +9.2 / -19.7 | Ba < 0, Bp +9 |

Bp follows the paper, but Ba falls below the hymen where the paper's stays above it. With one shared wall stiffness,
every change tried moves Ba and Bp nearly together, so the model meets either the paper's Bp or its Ba, not both
(README_2026-10-02.md).

## The earlier model: `L12_stabdamp_p3` (2026-09-24, superseded)

The most stable version as of 2026-09-24: it runs to full load (t = 1) with **no failed time steps**, in about 12 minutes on
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
  and runs in `claude_diag/runs/<name>/<name>.feb`. Tools that read results import the skill's `scripts/` from
  `~/.claude/skills/abaqus-febio-fea-pipeline/scripts` (copy the skill there first, as above). Many build scripts and tools default to earlier models in the
  chain, which are not included; the models above can serve as bases for new variants (most builders name their base in the
  script; `build_batch280.py` takes the source run on the command line). The plot files (`.xplt`) are too large for GitHub: run a model to regenerate its results.
- Every variant is written with a `<name>.feb.changes.txt` next to it that lists each change, and marks anything not in
  the Abaqus source as NOT IN SOURCE.

## Licence

No licence has been chosen yet, so all rights are reserved by the owner. Ask before reusing or redistributing.
