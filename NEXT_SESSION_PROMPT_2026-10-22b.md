Continue the PVP3D FEBio work in F:\Developing PVP3D Model for FEBio (hand-off 2026-10-22 ~15:50; the previous context filled
up). Work in claude_diag. The user is at the computer, steering the curve-end / strap work and viewing runs in FEBio Studio.

**Step 1: invoke the LOCAL skill `abaqus-febio-fea-pipeline`.** Then read claude_diag/README_2026-10-09.md section 8 from
"10-21 ~22:15: new session" to the end, NEXT_SESSION_PROMPT_2026-10-22.md (night + midday state), this file, and memory
(project-curve-selfpen-rules incl. its 10-22 correction, project-curve-pressure-fix, feedback-avw-damping-with-rest,
feedback-jitter-mass-damping, feedback-no-shared-cores, project-paper-sync).

**Runs: NO SHARED CORES**: `bash claude_diag/tools/queue_threads.sh 32 8 RUN` (4 x 8 threads); `tools/autostop_runs.sh 30 900 RUN`
for long runs (the queue's own watch is 300 min). Keep run paths <= ~255 chars (runs/<n>/<n>.feb.changes.txt), else shorten
names (cp0125 = curvepress0125, sc = strapcurve, sref = strap refined, h15 = hold to 1.5). Background processes survive Claude
restarts; background WAITERS do not (re-arm with tools/wait_any_end.sh 30 ... / wait_t.sh 30 T ...).

## The big correction of the day: how the curve is judged
- tools/curve_selfpen.py UNDER-REPORTS (node-in-hex depth to the nearest surface, capped by ~half the 1-2 mm curve wall; nodes
  that went right through are not counted). Use **tools/curve_cross.py RUN --times ...**: curve surface edges piercing curve
  faces, with depth (lower bound) and rest-position clusters. Every run crosses > 1 mm late: L350b (no curve pressure) ~3 mm
  distal; pressured runs 1-2 mm distal (rest z ~-27..-34) and ~1-1.4 mm at the mid-curve BEND (rest z ~-8..+2). The 10-20..10-22
  "curve <= 1 mm / clean" verdicts are withdrawn (README rows 14:45; paper addendum CORRECTION row).
- Curve vs other parts: a strap-crossing check (in README 15:30 row; code was a scratch script: curve_cross.surface + the strap
  domains' exterior faces) and /tmp/cross.py (curve nodes inside other parts) exist only as scratch; make a tool if needed.

## What the runs showed today (details: README section 8)
- Springs line L353 (stiff strap x10 tied direct + strap<->curve contact) + 1/8-1/4 curve pressure reaches Ba/Bp ~+10..+12 by
  t 0.88-0.92, then the AVW surges / folds back (L389, L391; L400 = L391 + curve_strap contact: same surge).
- AVW mass damping C 2000 (numerical aid; the user accepts it only with a settled hold, memory feedback-avw-damping-with-rest):
  L393 t 1 Ba/Bp +8.2/+9.3; in the hold the wall keeps descending (L394, L399c, L403: Ba/Bp ~+12/+20 at t ~1.06) and the LOWER
  curve folds into a U (rest z ~-12..-27, crossings 2-2.3 mm) -> crawl.
- L399c (curve nodes on all strap faces, curve_strap sr 0.25) / L403 (= L399c + strap refined 4x, build_batch404): curve self-
  crossings identical; but the finer strap halves the strap->curve penetration at the tip (0.5-0.9 vs 1.1-1.7 mm; the user: "kept it
  together better"). Curve material from rest z -6..-25 still passes through the strap from t ~1.0 (contact missed: sr 0.25 vs
  fast motion; some through the strap's back/edges, not in the contact). The user: it first pokes through between the strap and
  the curve where a curve element stretches a lot.
- USL x10 (L386): cuts the cervix's forward motion but raises it, less prolapse (the USL forces ~1-1.4 N a side either way: the
  source law has a 25 mm slack toe).

## Running at hand-off (15:46) and queued
| run | tests | state |
|---|---|---|
| `L403_springpatch_edge_fanstrap_strapmat10_edgedirect_sc_cp0125_MAXl_ant1_C2000_h15_cs_sref` | L399c + strap 4x finer | t 1.089 in the hold, folding (crossings 2.3 mm); stop when the user is done looking |
| `L402_springpatch_edge_fanstrap_strapmat10_edgedirect_strapcurve_cp0125_MAXl_ant1_flipcap` | L391 + strap flipped forward + strap-matrix cap over the curve end faces (build_batch403) | t 0.875 |
| `L401_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_cp05_avwp_MAXl_ant1` | L381 + extra AVW pressure 0.002 MPa fading 0.5 -> 0.8 (the user's 'push', corrected) | t 0.60, slow |
| `L388c_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_cp025_MAXl_ant1_hold15` | the user's HOLD on L384 (to t 1.5) + AVW-Para loft damping (L388 / L388b stalled on the loft jitter) | t 0.71 |
| queued: `L404_..._fillrim` (rim bevel, L391), `L409_..._infil20` (2 mm inner opening fillet, L391), `L410_..._C2000_h15_sref_if2` (L403 + 2 mm inner fillet), `L411_..._C2000_h15_sref_cr` (L403 + curve split 2x around the U below rest z -25, build_batch286) | | |
Judge each with curve_cross.py (and curve-vs-strap), POP-Q (pop_q_bump.py), node_speed.py; damped runs at a settled hold end.

## Built, not queued (the user looked)
L405 (lengthwise-crease fillet: "doesn't look correct"), L406 / L407 (opening fillet R 1.0 / 1.6), L408 (R 2.5 fillet + 1 mm
flipped strap: "does not look right"). Not built: rigid cap ABOVE THE AVW (L397 was wrongly above the cervix), the OUTWARD fillet
from the curve to the strap (needs added material, a gusset), the fillets on L402's geometry (the user: maybe later), a thicker
strap on the L403 setup.

## Tools added today
curve_cross.py (curve surface crossings), midline_bulge.py (wall midline outward by level), build_batch397 (rigid plane / AVW
traction / extra AVW pressure), 399 / 399c (strap faces in the curve self-contact / separate curve_strap contact), 402 (end-layer
merge: 0 mergeable pairs, tied), 403 (flip + cap), 404 (strap refine 4x), 405 (rim bevel / crease fillet), 406 (opening fillet with
harmonic morph, fold checks, --holdouter), 407 (flip + thicken strap, material x th0/th).

## Next ideas the user raised
Finer curve mesh where it stretches against the strap (L411); fillet spreading the stretch (L410); larger curve_strap search
radius (needs a thicker strap, sheet 0.49 mm); include the strap's back side / edges in the contact; the bend-region crossing
(rest z ~-8..+2) is untouched by all tip fixes.

## Before you finish
Update README section 8, run_notes (+ `py -3.10 tools/run_log.py` rebuild, ~20 min, no timeout), the paper addendum, memory
(project-pvp3d-handoff), and write the next NEXT_SESSION_PROMPT.
