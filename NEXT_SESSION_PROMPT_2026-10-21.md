Continue the PVP3D FEBio work in F:\Developing PVP3D Model for FEBio (hand-off 2026-10-21 morning, after the overnight session
on the CURVE COLLAPSE; the user was away and ask-before-runs was waived for that night only: from now on ASK before runs again,
unless the user waives it; 10-21 ~09:30 the user said "keep running files and learning" while away).

**Step 1: invoke the LOCAL skill `abaqus-febio-fea-pipeline`.** **Runs: NO SHARED CORES**: `bash claude_diag/tools/queue_threads.sh
32 8 RUN`; background waits `bash tools/wait_any_end.sh 30 RUN ...` with run_in_background, re-armed each time.
**Then read** claude_diag/README_2026-10-09.md section 8 from "2026-10-20 night (new session" to the end, this file, and memory
(feedback-curve-first-no-popq, project-curve-selfpen-rules, feedback-no-shared-cores, feedback-jitter-mass-damping, project-paper-sync).

## Standing instructions (unchanged)
- **Curve first, no POP-Q**: judge runs ONLY on the curve: self-penetration (tools/curve_selfpen.py RUN --times ...), tip kinks
  (tools/tip_kinks.py), the distal-curve span + PeB rotation (tools/curve_span.py, NEW 10-21), element volume in the bend, pace.
- Contacts "for later" are ONLY AVW-Para <-> PM arch and strap <-> vaginal wall.
- One change per variant first, then combine. Full run names in tables. Record every run in tools/run_notes.py + README section 8;
  paper decisions in claude_diag/PAPER_ADDENDUM_since_2026-10-19.md.

## What the night found (base L352 = springs line, original strap tied direct; collapses at t ~0.67-0.70 snap + 0.83)
(final numbers: README section 8 rows of 10-20 night / 10-21)

Deepest curve self-penetration L / R (mm) at t 0.70 and 0.83-0.85, the distal span (ratio to rest, right) and the end:

| run | the one change (vs base) | t 0.70 | t 0.83-0.85 | span 0.70 / 0.85 | end |
|---|---|---|---|---|---|
| `L352_springpatch_edge_fanstrap_edgedirect_MAXl_ant1` | base | 0.72 / 0.82 (snap at 0.67) | 0.89 / 0.91 (0.83) | 0.58 / 0.41 | AVW fold-back blow-up 0.883 |
| `L361_springpatch_edge_fanstrap_edgedirect_curvecontact_MAXl_ant1` | (1) +316 partly-tied distal faces in curve_all | 0.72 / 0.83 | 1.01 / 1.03 (0.838) | 0.58 / 0.42 | stopped 0.848: no help |
| `L362_springpatch_edge_fanstrap_edgedirect_prebend4_MAXl_ant1` | (2) distal canal pre-bent 4 mm +y (sin^2) | 0.28 / 0.35 (no snap) | 0.70 / 0.94 (0.79) | 0.56 / - | AVW blew up 0.796 (earlier) |
| `L360_springpatch_edge_fanstrap_edgedirect_sph010_MAXl_ant1` | (3) LA_sphincter_side_conn 0.25 -> 0.1 | 0.25 / 0.33 (no snap) | 0.48 / 1.01 (0.85) | 0.60 / 0.62 | stopped 0.876; PeB only +0.5 deg |
| `L357_springpatch_edge_fanstrap_edgedirect_strapcurve_MAXl_ant1` | strap <-> curve contact | 0.70 / 0.95 | 0.36 / 1.44 (0.86) | 0.58 / 0.38 | stopped 0.95 (crawl) |
| **`L359_springpatch_edge_fanstrap_edgedirect_curvepress_MAXl_ant1`** | **the curve pressure back on (ours, NOT IN SOURCE: the source has no curves)** | **0.23 / 0.23** | **0.32 / 0.44** | **1.20 / 1.14** | 0.84 / 0.75 at 0.90-0.91; distal posterior curve stalls the run at 0.925 |
| `L363_springpatch_edge_fanstrap_edgedirect_curvepress_sph010_MAXl_ant1` | L359 + (3) | 0.20 / 0.23 | 0.30 / 0.44 | 1.19 / 1.08 | curve gives way 0.887 (earlier than L359) |
| `L365_springpatch_edge_fanstrap_edgedirect_prebend4_curvepress_MAXl_ant1` | L359 + (2) | 0.25 / 0.23 | 0.28 / 0.55 | 1.25 / 1.21 | stopped 0.92: same late surge, a little earlier; dropped |
| `L366_..._curvepress_distalmass100_MAXl_ant1` | L359 + distal curve density x100 (jitter rule) | - | - | - | stopped 0.615: lag (PeB 8 vs 24.5 deg) |
| `L367_springpatch_edge_fanstrap_edgedirect_curvepress_distalmass10_MAXl_ant1` | L359 + distal curve density x10 | - | - | - | stopped 0.05 (not needed) |
| `L364_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress_MAXl_ant1` | full candidate L350b + curve pressure | 1.18 / 0.92 (no snap; tip-corner crush) | 0.77 / 0.60 (0.79) | 1.32 / - | stopped 0.7993: stuck on the right AVW-Para loft |
| `L350b_springpatch_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_MAXl_ant1` | the candidate (10-20) | 0.85 / 0.93 (snap at 0.69-0.70) | 0.76 / 1.08 | - | t 1, ~1 mm from 0.70 |

**Conclusions**
- **Putting back OUR curve pressure** (Load-curves: the 0.014 MPa wall pressure carried onto the curves' outer skin, 512 faces;
  added in L140, removed on the tube line in L254 / L258; the PVP3D source has NO curves and NO curve pressure).
  With it the distal curve is held stretched (span 1.1-1.25x rest instead of squeezed to 0.4x) and the collapse is gone through
  t ~0.9 (<= 0.55 mm). It is NOT a fidelity fix (corrected 10-21). Mechanism (found 10-21): the pressure squeezes the soft,
  nearly incompressible curves (Curve_soft = Vagina_AVW x0.6, ours), they lengthen and jack the apex ~22 mm ABOVE rest; at
  t ~0.88-0.93 that lift is released (cervix drops 42 mm: the "flies up, then slams down" the user saw) and the run crashes.
- Late (t > 0.88) the span still starts to shorten as the PeB rotates on, and the distal POSTERIOR curve (canal_wrap_hex_pvw /
  _pvw_fine) then stalls the run (L359 0.925, L363 0.887): the next problem.
- The pre-bend (2) only turns the 0.67 snap into a steady fold and brings the AVW fold-back forward; with the pressure it adds ~nothing.
- Softer sphincter-side connectors (3) barely free the PeB (+0.5 deg; +2 with pressure) and make the late curve failure earlier.
- More self-contact faces (1) and the strap <-> curve contact change nothing.
- The jitter mass x100 on the distal curve drags the PeB (the curve end is tied to it): rejected; x10 (L367) is the bracket.
- In the full candidate (L364) the pressure removes the U-fold snap, but ~1 mm remains as a LOCAL crush at the curve-tip corner where
  the direct-tied strap edge and the AVW-Para loft corner pull (rest (+-10..13, 12.9, -36.5), adjacent hexes); the run slows at
  t 0.80 on the right AVW-Para loft.


## The user's decisions 2026-10-21 morning
- "I think 359 fixed the issue, but then the prolapse was just too large": the curve is fixed by the curve pressure; POP-Q
  allowed again. L359 Ba / Bp / C -8.0 / -7.8 / -48.1 (t 0.85) -> +8.7 / +13.0 / -26.8 (0.92): a 17-21 mm late surge.
- Runs must go to t = 1 ("I still want to go to t=1. Just less impairments probably"): never stop at the old crash point.
- Test the curve pressure on less-impaired MAX-lite cases (all on L359, curve pressure kept):

| run | change vs L359 | status |
|---|---|---|
| `L368_springpatch_edge_fanstrap_edgedirect_curvepress_MAXl_ant1_la60` | levator 60 % (whole tissue) | running (07:24) |
| `L369_springpatch_edge_fanstrap_edgedirect_curvepress_MAXl_ant1_post0075` | P-arcus 0.05 -> 0.075 (local) | running (07:24) |
| `L370_springpatch_edge_fanstrap_edgedirect_curvepress_MAXl_ant1_la60_post0075` | both | running (09:21) |
| `L371_springpatch_edge_fanstrap_edgedirect_curvepress_MAXl_ant1_la70` | levator 70 % (whole tissue) | running (09:21) |

  Judge: curve <= ~1 mm to t 1 and no stall; POP-Q targets (Ba / Bp >= +10, apex motion >= 35 mm) at t 1; does the prolapse
  arrive gradually? Through t 0.6 L368 / L369 = L359 within 0.6 mm (the levator / P-arcus act late on the springs line).
- A new base only when one of these reaches t 1 with the targets; then BASES in run_log.py + default-bases memory + the three
  versions (memory feedback-three-base-versions). Ask the user first.
- Box: "the Box folder" = C:\Users\mgordon\Box\Research\Home Computer Results; "move a model" = COPY its .feb + .xplt
  (L359 copied 10-21; its .xplt was still uploading in Box Drive).

## Open questions for the user
1. Curve pressure (ours): fixes the collapse by jacking the apex up; then the slam. Next test proposed: L372 = L359 with the curves at the wall stiffness (Curve_soft x0.6 -> x1).
2. The candidate's ~1 mm tip-corner crush (L364): accept as within the ~1 mm rule, or treat it (e.g. the AVW-Para loft corner off
   the strap corner, or the pressure off the few tip faces under the strap / loft)?
3. The late (t > 0.88) distal posterior-curve stall: L367 (density x10) is the bracket; other levers: a finer distal mesh WITH the
   pressure (L286 folded without it), or a load hold.
4. Pre-bend and sphincter 0.1: dropped (L365 / L363 added nothing late).


## Update 2026-10-21 evening (the lift-and-slam; the user offline from ~18:30)
**Mechanism (README section 8 rows of 10-21):** our curve pressure (NOT in the source: it has no curves) squeezes the soft,
nearly incompressible curves (p / curve shear modulus 5.6); they LENGTHEN (centreline 1.33x at t 0.6) and the walls stretch
with them (wall edges 1.33-1.35x); the net curve-pressure force is small and slightly downward, so the cervix rise (+18-23 mm,
"the arch") is the lengthening pushing the cervix up between tied ends. At t ~0.88 the long curve gives way, the cervix drops
42 mm (the "slam"), the PVW rotates and its follower pressure turns downward (net wall load x5), the PVW twists (V/V0
0.48-2.17) and the distal posterior curve flies off. The lengthening is also what kept the distal curve stretched (no hairpin).
New tools: tools/cervix_lift.py (cervix u_y / u_z by t), tools/build_batch374 (CL / USL stop law), 375 (pressure on the distal
faces only).

| run | change vs L359 (= L352 + curve pressure 0.014) | result |
|---|---|---|
| L368 / L369 / L370 / L371 | levator 60 % / P-arcus 0.075 / both / levator 70 % | same arch and slam ~0.01 later; stopped |
| `L372_..._curvewall1_MAXl_ant1` | curves x1 of the wall (was x0.6) | **t 1 NORMAL**, no slam, curve <= 0.52 mm, but arch kept (+15 mm) and prolapse short (Ba / Bp -5.2 / -5.0) |
| L373 | USL x0.5 | MORE arch (+20.2 vs +18.4); stopped |
| L374 | CL / USL: x0.2 to 35 mm, then 3x the healthy end slope | arch peaks +14.6 (0.5) and recedes to +6.3 (0.82); crawled at 0.816; stopped |
| L375 | pressure on the distal faces only | no arch, but the hairpin is back (= L352); stopped 0.752 |
| `L376_..._curvepress05_MAXl_ant1` | curve pressure x0.5 (0.007) | **best balance**: arch +7.6-7.9, receding, no slam, curve <= 0.87 mm, Ba / Bp +8.6 / +16.2 at 0.90; then the anterior wall folds back (Aa -22 vs Ba -5 at 0.907); stopped |
| `L379_..._curvesrcwall_MAXl_ant1` | curves at the SOURCE vaginal wall law (x4.8; walls stay x0.35) | running: no arch to t 0.40 (-1.4 vs L359 +5.0); slow |
| `L380_..._post02` | P-arcus x0.2 | running: arch = L359 (+13.7 at 0.5) |
| `L381_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_MAXl_ant1` | full candidate L350b + curve pressure 0.007 | running (17:02); the user: the AVW-Para <-> PM arch interaction may stop L376's fold-back |
| `L382_fanstrap_parastrap_PM_band12_lip_tipanc_strapmat3_lipC500_edgedir_curvepress05_csw_MAXl_ant1` | L381 + curves at the source wall law | running (18:26) |

**When the user reconnects:** report L379 / L380 / L381 / L382 with tools/cervix_lift.py (arch), curve_selfpen.py (curve),
pop_q_bump.py (Aa vs Ba for the fold-back, targets Ba / Bp >= +10), element_volume_report _PickedSet64 (PVW twist), and
whether each reached t 1. Ask before new runs. Other user notes 10-21: CL / USL are CONN3D2 force-elongation tables (three
tables, no area; figures/cl_usl_spring_law.png); the user reads ~1 MPa for US and vaginal wall (ours: walls 12.5 kPa at small
strain, curves 7.5; source wall 50 kPa) and is looking into it before any modulus-based law; Box folder rule (memory).
- The user (10-21 ~22:20): in the coming 8-10 h run, do one HOLD run: the best candidate (L384 if good, else L381) with time_steps
  100 -> 150 (full load held t 1-1.5; settle_damping C 20 switches on at t 1.0-1.05; it is OFF during loading). Report Ba / Bp / C
  and the curve from t 1 to 1.5.
