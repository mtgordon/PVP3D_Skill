Continue the PVP3D FEBio work in F:\Developing PVP3D Model for FEBio (hand-off 2026-10-23 morning, written while the user was
away overnight). Work in claude_diag. Step 1: invoke the LOCAL skill `abaqus-febio-fea-pipeline`, then read
claude_diag/README_2026-10-09.md section 8 from "10-22 ~16:00" to the end, this file, and memory (project-canal-width-source,
feedback-finer-strap-mesh, project-curve-selfpen-rules, feedback-avw-damping-with-rest, feedback-no-shared-cores).

## The result of the night: the canal width
- The user's idea (10-22 evening): our walls + curves (+ cervix) were 8-28 % WIDER than the source AVW / PVW / cervix
  (PVP3DModel_job.inp), so the wide canal is dragged through the hiatus and folds. tools/build_batch422.py narrows the canal
  (AVW, PVW, cervix, apex_contour, curves) in x by level to the source width (shrink only) and spreads PM_conn's 26 springs IN
  ORDER over the whole distal edge (27 AVW nodes + the curve corner each side; no crossing). Picture claude_diag/canal_width_L422.png,
  claude_diag/pmconn_L423.png.
- **L423b** (= L403's smaller fine strap + the width; L423 had a duplicate spring_patch_ties block, fixed in build_batch288):
  the first run to NORMAL TERMINATION at t 1.5, but the hold collapse is only DELAYED to t ~1.45-1.48 (the user: "poking through for a while and then collapsing around 1.48"; settled to 1.4 with a steady ~1.1 mm distal-curve crossing at rest z -31, then it jumps to 2.4 mm and the walls follow). Curve self-crossing ~0.9-1.1 mm to t 1.2
  (was 1.9 mm and growing in L413b / L403), walls never through the PM arch. BUT the damped wall keeps creeping through the
  hold (node speed 0.22 mm per 0.001 t at 1.5), Ba / Bp -2.7 / -3.9 (t 1) -> +0.1 / -0.1 (t 1.5), and the crossings grow
  slowly to 1.9 mm by t 1.5 (left distal curve, rest (-6.6, 10.5, -30)).
- The narrow canal prolapses much less: MAXl_ant1 Ba / Bp -2.7 / -3.9 at t 1 (L413b +7.7 / +9.4); **MAX (L428) only -2.5 /
  -3.5 at t 0.96**: the range targets (Ba / Bp >= +10) are NOT met any more. Healthy (L423h) -18.8 / -21.5 (fine).

## Contact tests (one change each)
| run | change | verdict |
|---|---|---|
| L416 (on L413b) | curve_strap penalty 5 -> 25 | halves curve->strap to t 0.9, hold collapse unchanged |
| L417 | curve self-contact penalty 20 | crawls from the start: unusable |
| L418 | self-contact offset 0.2 | stopped early (slot) |
| L419 / L420 | curve_strap sr 0.4 / all strap faces in curve_strap | stopped early (L416 showed the contact side cannot stop the fold) |
| L414 | self-contact only on outer + distal faces | = L413b |
| L421 (L413b) / L426 (L423b) | canal <-> PM arch contact (build_batch421) penalty 1 sr 1 / penalty 5 sr 2 | weak one holds nothing; strong one folds the left distal curve (1.83 mm): NOT usable as built (coarse 70-facet one-sided PM) |
| L427 (L423b) | curve_strap penalty 25 | strap poking halved (0.5-0.8 mm), but left distal curve self-crossing worse (1.84 mm at t 1) |
| L411 / L415 / L413 | curve mesh refinement (around the U / along the canal vaginal part / whole curve) | no help; whole-curve split crawls (cervix / apex ties) |
| L409 / L410 | 2 mm inner opening fillet | crawl from the start |
Tools added: strap_cross.py (strap vs curve crossings), build_batch412 (along-canal curve split, --zmax), 413 (strap widened:
the user's red area), 414 (self-contact face classes), 416 (contact settings), 420 (all strap faces in curve_strap), 421
(canal <-> PM contact), 422 (canal width to source + PM_conn in order), 424 (support case switch from ANY case; healthy / MAX).
PM crossing check = scratch pmcross.py (canal surface edges through PM_Plane); make a tool if needed.

## Still running at hand-off (check them first)
L427 (L423b + curve_strap 25, MAXl_ant1, hold), L428 (L423b MAX), L429 (L427 MAX), L430 (L427 healthy).

## Open for the user
1. The narrow canal no longer reaches the MAX range targets: accept, or find what else (the source had no curves; maybe the
   curves now carry load the source walls did not?).
2. The creep in the hold (damping C 2000 slows the settle; memory feedback-avw-damping-with-rest wants a settled end): longer
   hold or less damping on L423b.
3. The left distal curve cluster (rest (-6.6, 10.7, -32)) is the remaining self-crossing in L423b at 1.5, L426, L427.
4. The solid AVW-Para back + its PM contact (the user's "next") was NOT ported: the L342-L344 chain is multi-step and crawled
   on the springs line; the canal <-> PM contact needs a finer PM or another form.
5. L404 (rim bevel) was taken out of the queue (coarse strap).

## Before you finish
README section 8, run_notes (+ `py -3.10 tools/run_log.py` rebuild, ~20 min), the paper addendum, memory, next NEXT_SESSION_PROMPT.
