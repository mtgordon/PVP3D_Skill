Continue the PVP3D FEBio work in F:\Developing PVP3D Model for FEBio (hand-off 2026-10-21 ~22:15; the previous context ran out
of room). **The user is away 8-10 h and said: run freely** (defaults agreed 10-21 ~22:10, below). Work in claude_diag.

**Step 1: invoke the LOCAL skill `abaqus-febio-fea-pipeline`.** Then read claude_diag/README_2026-10-09.md section 8 from
"2026-10-20 night (new session" to the end (the whole 10-20 night + 10-21 story), NEXT_SESSION_PROMPT_2026-10-21.md (its tables
and the "Update 2026-10-21 evening" section), this file, and memory (project-curve-pressure-fix, project-max-lite,
feedback-no-shared-cores, feedback-jitter-mass-damping, project-paper-sync, reference-box-folder, feedback-curve-first-no-popq).

**Runs: NO SHARED CORES**: `bash claude_diag/tools/queue_threads.sh 32 8 RUN` (4 runs x 8 threads); `tools/autostop_runs.sh 30
MAXMIN RUN` watches stalls only for MAXMIN minutes: L383 / L384 / L385 have 300-min watches started 20:36-21:24 -> RE-ARM them
(e.g. `nohup bash tools/autostop_runs.sh 30 900 RUN &`). Waits: `bash tools/wait_any_end.sh 30 RUN ...` / `tools/wait_t.sh 30 T RUN`
with run_in_background, re-armed each time (keep one pending while unattended). Long run names: the .feb.changes.txt can fail
(Windows 260-char path); then the README row is the change note.

## The user's go for this 8-10 h stretch (defaults agreed 10-21 ~22:10)
1. Start / stop / combine runs freely within 32 threads; stop only runs that are stuck or already answered; log everything
   (tools/run_notes.py, README section 8, claude_diag/PAPER_ADDENDUM_since_2026-10-19.md).
2. Do one HOLD run: the best candidate (L384 if its prolapse is good and the curve < 1 mm, else L381) with `time_steps` 100 -> 150
   (full load held t 1 -> 1.5; settle_damping C 20 switches on at t 1.0-1.05, it is OFF during loading). Name `..._hold15`.
   Report Ba / Bp / C and the curve from t 1 to 1.5.
3. Success: curve self-penetration <= ~1 mm to t 1; no slam; no anterior fold-back (Aa ~ Ba); no PVW twist (V/V0 within ~0.6-1.6);
   a small arch (cervix within a few mm of rest; the user: the L379 arch, +7 mm, is NOT small); Ba / Bp >= +10, apex motion >= 35 mm.
   POP-Q is judged again (the curve is mostly fixed).
4. Combine what helps (e.g. L386's stiff USL with quarter pressure / MAX), one change per step, carry the best to t 1.
5. NOT approved: the "stiff curves only on the distal curve, soft above" idea is MINE (Claude's), not the user's: do NOT build it
   unless the user says yes.

## The physics found 10-21 (details in README section 8 rows of 10-21)
- The curves and their pressure are OURS (the PVP3D source has no curves): Load-curves = the 0.014 MPa wall pressure carried onto
  the curves' outer skin; Curve_soft = Vagina_AVW x0.6; the walls are x0.35 of the source vaginal wall law.
- With the curve pressure the soft, nearly incompressible curves are squeezed and LENGTHEN (1.3-1.5x); the walls stretch with them;
  the cervix rises / moves forward ("the arch"). Net curve-pressure force is small. At t ~0.88 the arch can release in a "slam"
  (L359). The lengthening is also what kept the distal hairpin away (pressure on the distal faces only, L375: no arch, hairpin back).
- **The curves act as RAILS**: with pressure (graded with the pressure) or stiffness, the canal moves out as a rigid tube (15-19 mm
  at every level at t 0.76) instead of the mid-canal bulging out 33-37 mm as in L352 (no curve pressure): the extra length goes
  into the arch instead of the prolapse. Not the introitus (moves out ~13-15 mm in all runs), not the strap (pulls hardest in L352),
  not the AVW-Para <-> PM contact (<= 1.7 N) (tools/lower_canal_hold.py).
- In the user's FEBio Studio view (+Z screen-left, +Y screen up-right) the visible "arch going up" is mostly the cervix's FORWARD
  motion (+y ~37-44 mm, +z ~10); the USL (fan USL-L/R_conn_2: 8 springs to one sacral hub) pull straight back against it. USL x0.5
  (L373) cut the forward motion 8 mm but added 2 mm of lift.
- At t 1 the runs are nearly settled (Amp-1 is a SMOOTH STEP): L381 creeping ~0.8 mm per 0.05 t, L379 at rest, L372 still descending.

## Runs at hand-off (22:11) and the key earlier ones
| run | change | state |
|---|---|---|
| `L381_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_MAXl_ant1` | full candidate L350b + curve pressure 0.007 (half) | **t 1 NORMAL (2 h 59 min)**: curve <= 0.52 mm, no slam / fold-back / twist, tip straight; arch +10 -> +7.9 mm; Ba / Bp / C -9.8 / -9.5 / -37.7 (short) |
| `L379_springpatch_edge_fanstrap_edgedirect_curvepress_curvesrcwall_MAXl_ant1` | L359 + curves at the source wall law | t 1 NORMAL: curve <= 0.16 mm; arch +7 (the user: not small); Ba / Bp -12.7 / -13.8 |
| `L382_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_csw_MAXl_ant1` | L381 + curves at the source wall law | running t 0.91: arch +4.0 at 0.76, curve <= 0.15, Ba -16 at 0.77 |
| `L383_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_MAX` | L381 at full MAX (AVW-Para fibre E x0.5) | running t 0.70: arch = L381 so far |
| `L384_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress025_MAXl_ant1` | L381 with curve pressure 0.0035 (quarter) | running t 0.55: cervix -2.6 at 0.4 (L381 +0.1) |
| `L385_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_csw_MAX` | L382 at full MAX | running t 0.27 |
| `L386_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_usl2_MAXl_ant1` | L381 + USL scale 0.2 -> 2 (2x healthy; the user's idea: pin the cervix so the length goes out) | QUEUED (queue_threads process waiting for 8 threads) |
| L376 (springs, half pressure, soft curves) | | best springs-line balance: arch +7.6, Ba / Bp +8.6 / +16.2 at 0.90, then anterior fold-back; stopped 0.907 |

Tools added 10-20 / 10-21: curve_span.py, cervix_lift.py, lower_canal_hold.py, build_batch360 (discrete-material scale),
361 (more curve-contact faces), 362 (pre-bend), 364 (add the curve pressure), 374 (CL / USL stop law), 375 (distal-only pressure);
build_batch356 has --df (density). Figures: claude_diag/curve_end_view/L352_vs_L359_curve_side.png, L362_prebend_rest.png,
claude_diag/figures/cl_usl_spring_law.png, cl_usl_law_L374.png. The user reads ~1 MPa for USL / vaginal wall in the literature
(ours far softer; README row 'vaginal wall modulus vs ~1 MPa'): looking into it, no modulus-based law yet.

## Before you finish
Write NEXT_SESSION_PROMPT_2026-10-22.md; update README section 8, run_notes (+ `py -3.10 tools/run_log.py` rebuild if the
spreadsheet is due), the paper addendum, and memory (project-pvp3d-handoff, project-curve-pressure-fix).
