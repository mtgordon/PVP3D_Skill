Continue the PVP3D FEBio work in F:\Developing PVP3D Model for FEBio (hand-off 2026-10-22 morning, after the 10-21 night
stretch the user ran unattended, 22:10 -> ~07:00). Work in claude_diag.

**Step 1: invoke the LOCAL skill `abaqus-febio-fea-pipeline`.** Then read claude_diag/README_2026-10-09.md section 8 from
"10-21 ~22:15: new session" to the end (the night's rows), NEXT_SESSION_PROMPT_2026-10-21b.md (the physics found 10-21: the
curve pressure is OURS, the pressurised curves act as RAILS), this file, and memory (project-curve-pressure-fix,
project-max-lite, feedback-no-shared-cores, feedback-jitter-mass-damping, project-paper-sync, feedback-curve-first-no-popq).

**Runs: NO SHARED CORES**: `bash claude_diag/tools/queue_threads.sh 32 8 RUN`; `tools/autostop_runs.sh 30 900 RUN` (the queue's
own watch is only 300 min; re-arm long runs). Waits: `tools/wait_any_end.sh 30 ...` / `tools/wait_t.sh 30 T ...`.
**Path length**: runs/<n>/<n>.feb.changes.txt must stay <= ~255 chars or FEBio cannot write the .xplt (L388 died at start,
261 chars) / the builders cannot write the change note (L392: README row only). Shorten names (cp025 = curvepress025,
sc = strapcurve, h15 = hold15).

## What the night found (details: README section 8 rows of 10-21 night / 10-22)
1. **The candidate line (full candidate L350b: strap x3 tied direct, AVW-Para fibre lofts, AVW-Para <-> PM contact) is
   capped**: even with NO curve pressure it stops at Ba +1.0 at t 1 (L350b, curve ~1 mm). Quarter curve pressure (L384) is its
   cleanest run (t 1: no arch, curve <= 0.78 mm, Ba / Bp -6.1 / -6.3). Full MAX adds ~1 mm (L383), stiff curves remove the arch
   but hold the canal (L382 -14.8), USL x10 (the user's idea, L386) cuts the cervix's forward motion 17 mm but RAISES it and
   gives 6 mm less prolapse, the strap <-> curve contact does nothing there (L390 = L384 within 0.1 mm), 1/8 pressure (L387)
   +5 mm by t 0.73 then stuck on a left AVW-Para loft jitter.
2. **The springs line L353 (stiff strap x10 tied direct + strap <-> curve contact; Ba +14 at t 1 but curve 0.9-1.4 mm) + a
   little curve pressure reaches the targets with a clean curve**: L389 (quarter) Ba / Bp +10.0 / +15.2 at t 0.92, curve <= 0.64;
   L391 (1/8) +12.2 / +15.9 at 0.88. **Then the anterior wall SURGES and folds back** within 0.01 t (AVW node speed x10-20,
   Aa -18 / -44), the run crawls and dies. Less pressure = faster prolapse = earlier surge.
3. **AVW mass damping C 200 -> 2000 (NUMERICAL AID) removes the surge** but lags the wall: L392 (quarter) t 1 Ba / Bp +3.6 /
   +5.5; **L393 (1/8) t 1: Ba / Bp / C +8.2 / +9.3 / -28.4, curve <= 0.69 mm, no arch (apex ~41 mm), no fold-back, PVW V/V0
   0.80-1.30, the mid canal BULGES (not a tube)**; both still creeping at t 1 (p95 ~0.1-0.12 mm per 0.001 t).
4. New tool tools/midline_bulge.py RUN --t ...: wall midline outward by level (tube vs bulge; reproduces the 10-21 rows).

## Runs at hand-off (check them first)
| run | change | state at ~06:50 |
|---|---|---|
| `L388b_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_cp025_MAXl_ant1_hold15` | **the user's hold run**: L384 + time_steps 100 -> 150 (held t 1 -> 1.5); rerun of L388 (stuck at 0.781 on a non-deterministic AVW-Para-L loft jitter) | running t ~0.5 (started 05:19; ~6 h to 1.5) |
| `L394_springpatch_edge_fanstrap_strapmat10_edgedirect_strapcurve_cp0125_MAXl_ant1_avwC2000_h15` | L393 held to t 1.5 (diagnostic: damping lag or equilibrium?) | running (06:23; ~3 h) |
| `L395_springpatch_edge_fanstrap_strapmat10_edgedirect_strapcurve_cp025_MAXl_ant1_avwC600` | L389 with AVW damping C 600 | running (06:40) |
| `L396_springpatch_edge_fanstrap_strapmat10_edgedirect_strapcurve_cp0125_MAXl_ant1_avwC600` | L393 with AVW damping C 600 (= L391 + C 600) | running (06:47) |
Finished / stopped tonight: L382 (t 1), L383 (error 0.993), L384 (t 1), L385 / L386 / L387 / L390 (stopped), L389 / L391
(surge), L392 / L393 (t 1).

## Next (to decide with the user)
- **Read L394 at t 1 -> 1.5**: if the held L393 settles at Ba / Bp >= +10 with no fold-back, L393 (+ its hold) is the
  candidate; if it surges once the load is held, the damping only delays the fold-back.
- L395 / L396: the least damping that avoids the surge.
- The user's hold on L384 (L388b) answers how far the candidate settles after t 1 (L384 was creeping ~0.6 mm per 0.05 t).
- **Ask the user**: is AVW mass damping C 600-2000 acceptable as a numerical aid in the reported model (it changes the
  late dynamics, not the materials), or should the surge be fixed physically (the candidate's AVW-Para <-> PM contact stopped
  the fold-back in L381, but that line is capped at Ba ~+1)?
- NOT approved, do not build: the 'stiff distal curve only' idea (Claude's).

## Before you finish
Write the next NEXT_SESSION_PROMPT; update README section 8, run_notes (+ `py -3.10 tools/run_log.py` rebuild), the paper
addendum, and memory (project-pvp3d-handoff, project-curve-pressure-fix).

## Update 2026-10-22 midday (the user back, steering the curve-end / strap work)
- Damping alone does not fix the springs line: L395 / L396 (AVW C 600) and the held L394 (C 2000) all run past Ba +10 and
  surge / fold back; the user accepts AVW damping only with a settled hold (memory feedback-avw-damping-with-rest).
- L394's curve went > 1 mm at the distal end corner (t 1.03): the strap is tied by its base along the end-face rim and was
  extruded BACK into the curve (overlaps the end hexes); the corner hexes are a column of slender, heavily tied elements
  (0.45 x 1.5 x 1.0 mm; junction / patch / end ties) that swell 2.5x / crush 0.3x; the end-face rim folds 71-94 -> 143-168 deg.
- Strap / curve-end runs (all NOT IN SOURCE): L399c (curve nodes on all strap faces, sr 0.25, damped + held) and L400 (same
  on L391) running; L403 (L399c + strap refined 4x; tools/build_batch404.py) running; queued: L402 (strap flipped forward +
  a strap-matrix cap tied over the end faces, own cap contact sr 0.3; build_batch403.py), L388c (the user's hold on L384 +
  AVW-Para loft damping C 200; L388 / L388b stalled on the AVW-Para loft jitter), L404 / L405 (fillets on L391: end-face rim
  / lengthwise creases; build_batch405.py; figure claude_diag/figures/curve_end_fillet_L404_L405.png).
- Stopped / not built: L397 (rigid cap plane placed above the cervix, the user meant above the AVW: rebuild as a sheet shaped
  over the AVW, or flat just above the AVW top), L398 (-z push; wrong direction), L401 = the corrected push (extra AVW
  pressure 0.002 MPa fading 0.5 -> 0.8) running; the corner re-mesh (build_batch402.py, 0 mergeable pairs: tied nodes).
- REMEMBER (the user): do the fillets on L402's geometry later; the thicker strap (~1.0 mm, fibre E and matrix c1 x0.49 to keep
  the membrane stiffness) on L402's flipped geometry is the next strap idea if L402 still lets nodes through.
- 14:20 update: L399c reached t 1.035 in its hold (Ba / Bp +11.3 / +16.1; the strap contact kept the curve <= 0.87 mm), then the
  AVW sped up and it failed: the curve is no longer what breaks it, the hold's runaway descent is. L403 (L399c + 4x strap) running.
  Built for the user to look at, NOT queued: L406 / L407 (opening fillet R 1.0 / 1.6) and L408 (opening fillet R 2.5 + strap flipped
  forward 1.0 mm, material x0.49; build_batch406 / 407). Queued: L402, L388c, L404. The OUTWARD fillet (curve -> strap) still to
  design (needs added material, e.g. a gusset). Background waiters must be re-armed after any restart.
