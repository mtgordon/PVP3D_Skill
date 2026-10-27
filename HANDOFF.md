# Hand-over log (newest first)

New here? Read `notes/PROJECT_BRIEF.md` first.

Each entry: date, who, what ran (full run names) and what it showed, what is still running and where, what is next, and who
has the turn. Only the person who has the turn starts runs.

## 2026-10-26, the colleague (long session 10-24 -> 10-26 noon; all runs on the colleague's PC)
- **Turn:** passes to Mark once the colleague pushes this (her call). Nothing running. Detailed hand-off: `NEXT_SESSION_PROMPT_2026-10-26.md`; journal: README section 8.
- **Curve corner self-crossing fixed.** The curve hexes in the fold zone (rest z -45..-27) at the curve law x1.5
  (build_batchC003). **Minimal fix = `C017_C004_band15_tip`**: corner 0 edges > 1 mm (was 1.42), strap 0, PM arch 0,
  swelling 2.75x (was 4x). Checked in all three cases: MAXl_ant1 C017 Ba/Bp -10.1/-12.0, MAX `C021_L470_LA91_band15tip_mp`
  -9.3/-11.3, healthy `C020_C002_band15_tip_mp` -21.1/-24.1 (unprolapsed). Band x1.25-x2, tip-only to rest z +10: all work.
- **Levator:** the colleague's L91 law + fixed arcus + LA contact settings (build_batchC001; C004 = L465 + levator) removes
  L465's PM-arch crossing and halves the strap crossing; nothing through the levator (max 8.4 mm).
- **Price:** less descent than C004 (-6.2/-7.2), but that extra descent was the corner collapse: softening the rest of the curve
  (C019) or the strap (C014) gives nothing back. MAX barely drops more than MAXl_ant1 (C012, C021), as in Mark's L470.
- **Left open:** a sub-mm mid-canal crossing = the curve skin through the curve's OWN end faces tied into the walls
  (curve_cross_tri_surfaces.py; the curve crosses no wall-only face, curve_wall_cross.py). End-face contacts do not act (C018,
  C022); a curve-vs-walls_LA_primary contact (C023) chatters because that surface contains the curve skin. Optional:
  `C024_C017_mid2` (+ curve x2 at rest z -12..+10) cuts it (> 0.5 mm edges 16 -> 4, deepest 0.90 -> 0.74) for ~0.5 mm less
  descent. The range-of-motion targets (+10 mm) remain far off (strap-layout limit, as Mark found).
- New tools: build_batchC001 / C003 / C004 / C008 / C018 / C023, curve_cross_diag.py, curve_cross_tri_surfaces.py,
  curve_wall_cross.py, mrb_fibres.py (Slicer .mrb fibre mapping, ready; the Michigan file has not arrived).
- Output: C-runs from C004 on use must-point output (~145 MB .xplt instead of 2-2.5 GB; plot_stride is ignored by FEBio).
- Bases in the repo (git add -f): C001-C024 .feb + .changes.txt (L26 is needed by pop_q_bump.py; symlink L149 if missing).

## 2026-10-24, Mark (overnight session): turn passes to the colleague
- **Turn:** the colleague. **Your focus (Mark): keep the curve corner from penetrating itself.**
- **New base: `L465_L431_ogmid`** = L431 (strap x1/3, source-width canal) with the AVW / PVW / curves in a "middle" Ogden law
  (tension = source, compression 1.6x softer). NORMAL TERMINATION t 1.5, settled; MAXl_ant1 Ba/Bp -2.7 / -4.0; healthy (L469)
  stays unprolapsed (all MAXl_ant1 = MAX-lite unless named; the full MAX case L470 ~-4.3 / -5.1 at t 0.95). Problem: the curve still gets caught on itself (~1.4 mm at rest z -31 left / -21 right; L431 ~1.24), the
  curves cross the PM arch by 0.75 mm, the distal curve elements swell to 4x. The softer-compression wall did NOT fix the corner.
- Also finished: L459, L466, L467, L468 (49 springs from the distal-face edge to the PM arch instead of the strap): clean but the
  AVW "pushes back instead of out" (Mark) - little prolapse. L473 (a temporary push on the AVW) built, not run.
- AVW-Para connectors: use the hybrid (patches where the curve is coarse, single attachment-line nodes in the refined distal zone;
  build_batch456 / L466); not the single nodes everywhere (crush) nor the old patches everywhere (stall / swell).
- Ideas for the corner (details in `NEXT_SESSION_PROMPT_2026-10-24.md`): a local stiffer band in the folding curve elements;
  pre-bend / a wider bend radius; refine around the fold; self-contact on the outer skin + distal faces at penalty 5; spread the
  strap's tip attachment; an explicit-dynamics check.
- Still running on Mark's machine at hand-off: L469 (L465 healthy), L470 (L465 MAX), L471 (L465 + softer Ogden), L472 (L466 MAX).
  Their .feb files are in the repo; Mark will add their results.

## 2026-10-23, Mark: the repo becomes the hand-over repo
- **Turn:** Mark (runs still going on Mark's machine).
- **Bases in this repo** (`claude_diag/runs/`): L403 (fine strap), L423b (canal at the source width), L431 (+ strap x1/3:
  first clean t 1.5 finish, old curve material), L439 (fibre curves + strap-side contacts; stalled t 0.64), L440g (+ fine
  distal curve + strap merged onto it), L444 (L440g line + 4-fibre AVW / PVW wall, source k), L445 (L444 + PM_conn ends on
  the AVW only), L446 (L444 + strap x1/3), L447 (L444 + stiffer wall matrix c1 0.0005, E 0.011), L448 (L446 + PM_conn on every AVW distal-row node).
- **Running on Mark's machine:** L444, L445, L447, L448 (L440g stopped at t 0.53: looked bad; L446 stopped for L448: it kept
  the old PM_conn layout). Results go into README section 8 and `run_log_notes.csv`.
- **Learned today:** a 2-fibre (0 / 90 deg) wall with a low bulk modulus blows up (L442c); use 4 in-plane fibres + the source
  k (`build_batch444.py`). The curves' fibre material is still 2-fibre: probably the same weakness (next to test).
- **Next:** judge L444 / L445; then the softer strap (x1/3, as L431) on top; 4-fibre curves.
- Details: `NEXT_SESSION_PROMPT_2026-10-23b.md` and the next one.
