Continue the PVP3D FEBio work (hand-off 2026-10-26 noon, after the colleague's long session 10-24 -> 10-26; all C-runs ran
on the colleague's PC through FEBio Studio). Work in claude_diag. The person whose turn it is steers; ASK before new runs.

**Step 1:** `git pull`; invoke the project skill `abaqus-febio-fea-pipeline`; read the top entry of HANDOFF.md, this file,
claude_diag/README_2026-10-09.md section 8 from "10-24" to the end, and notes/MEMORY_NOTES.md.

## What the colleague's session did (Mark's email of 10-24: levator in all models, fix the curve corner, Slicer fibres)
1. **Levator package** (`build_batchC001.py SRC CASE NEW`; CASE = the source run's support case from build_batch244.CASES):
   the colleague's L91 Yeoh levator (c1 0.0469 / c2 0.0657 / k 4 = 2x LA100, at the case percentage; NOT FAITHFUL), the 62
   ATLA / P-arcus chain nodes 23677-23738 fixed (NOT IN SOURCE), PVW_LA / walls_LA penalty 0.5 -> 2, search_tol 0.01 -> 0.1
   (NOT IN SOURCE). On L465 (C001 / C004): removes the PM-arch crossing, halves the strap crossing, levator max 8.4 mm.
   build_batch424 asserts the old LA law: switch cases BEFORE applying C001.
2. **Plot output** (`build_batchC004.mp_output`): PLOT_MUST_POINTS + a constant dtmax load curve with a point every 0.01
   -> ~151 states, ~145 MB .xplt (L465-line runs wrote 2-2.5 GB; FEBio ignores plot_stride).
3. **Curve corner fix** (`build_batchC003.py SRC NEW --zlo --zhi --scale --tag`): curve hexes whose rest centroid is in
   [zlo, zhi) go to twin domains `<parent>_<tag>` with `Curve_<tag>` = Curve_ogmid c1 / c2 / k x scale (same exponents) and the
   parent's damping. **Minimal fix = x1.5 on [-45, -27)** (C017). The fold sat at the curve tip where the strap attaches;
   a band only in the middle of the fold zone (C005) moved the fold to its lower edge; softening the strap at the tip (C008)
   did nothing.

## Results at t 1.5 (MAXl_ant1 unless named)
| run | change | corner | mid-canal | strap | Ba / Bp |
|---|---|---|---|---|---|
| L465 (Mark) | base | 1.41 | - | 1.08 / 0.81 | -2.7 / -4.0 |
| C004 | L465 + levator | 1.42 | 0.84 | 0.40 / 0.65 | -6.2 / -7.2 |
| C010 / C011 / C016 | band x1.5 [-45,-19) / x1.5 [-45,+10) / x1.25 [-45,+10) | < 1 | 0.75-0.84 | 0 / 0.47 | -9.1 .. -10.4 / -11.0 .. -12.1 |
| **C017** | **band x1.5 [-45,-27)** | **< 1** | 0.85 | **0** | **-10.1 / -12.0** |
| **C021** (MAX) | L470 + levator + C017's band | < 1 | 0.85 | 0.41 / 0.60 | -9.3 / -11.3 |
| **C020** (healthy) | C002 + C017's band | < 1 | 0.99 | 0 | -21.1 / -24.1 |
| C019 | C017 + rest of curve x0.75 | < 1 | 0.90 | 0 | -10.1 / -12.0 (no descent back) |
| C024 | C017 + curve x2 [-12,+10) | < 1 | 0.74 (4 edges > 0.5 vs 16) | 0 | -9.4 / -11.7 |
("< 1": no curve-crossing edge deeper than 1 mm.) Stopped: C007 (x4, crawled), C008, C014 / C015 (strap softer = no change),
C018 / C022 (end-face contacts, never act), C023 (see below).

## Conclusions
- The corner fix works in all three cases, keeps the healthy model up and removes the strap crossing.
- The ~4 mm less descent than C004 is the corner collapse being removed, not a cost one can buy back (C019, C014).
- MAX barely drops more than MAXl_ant1 (C012, C021): the range (+10 mm) is a separate, strap-layout problem (Mark's L470).
- The remaining sub-mm mid-canal crossing is the curve skin through the curve's OWN end faces tied into the walls
  (`curve_cross_tri_surfaces.py`); the curve crosses no wall-only face (`curve_wall_cross.py`, max 0.03 mm). NOTE:
  walls_LA_primary contains the whole curve skin (all 1910 curve_all nodes), so a contact against it (C023) puts the curve
  against itself at zero gap and chatters. C024's second band cuts the crossing; adopting it is optional.

## Suggested next
1. Decide whether C017's band (or C017 + C024's band) becomes part of the base for the strap-layout work on the range.
2. The range target: MAX is still Ba / Bp ~ -9 / -11 vs the target ~ +10 (Mark's strap-layout question).
3. Slicer fibres: `claude_diag/tools/mrb_fibres.py` (inspect / map: VTK, STL, mrk.json, NRRD readers, ICP to the model's
   levator, per-element mat_axis ElementData) is ready and tested on a synthetic scene; the Michigan .mrb has not arrived.

## Practical notes
- pop_q_bump.py reads L26_springs_la3_pm (not in the repo): symlink/copy L149_seamspr_pm_outer as L26 if missing.
- The C-runs were started through FEBio Studio (Edit XML, "Save model before running" unticked), 2 at a time.
