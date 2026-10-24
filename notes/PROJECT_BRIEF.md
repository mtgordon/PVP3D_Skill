# PVP3D in FEBio: what we are doing, why, and how (brief, 2026-10-23)

## What
PVP3D is a 3D pelvic-floor model built in Abaqus (`PVP3DModel_job.inp`): the vaginal canal (anterior wall AVW, posterior wall PVW,
cervix), the levator ani (LA), the perineal body (PeB) and perineal membrane (PM), and the ligaments / attachments as 291 `CONN3D2`
connectors (cardinal CL, uterosacral USL, paravaginal AVW-Para, arcus, PM connectors). It is loaded with ramped abdominal
pressure, and impairing the supports makes the walls prolapse. The aim is to run it in FEBio 4.13 (open source, so the model and
the pipeline can be shared and used for parameter studies / ML later), and to show it reproduces the paper's range of motion
(Luo et al. 2015, J Biomech: Ba / Bp >= +10 mm below the hymen and >= 35 mm apex motion in the impaired cases, while the healthy
model stays unprolapsed), with POP-Q points measured as Bump et al. 1996 define them.

## Why it is hard
- The source runs `*Dynamic, Explicit`: no equilibrium iterations, so near-mechanisms, loose seams and buckling parts just ride
  through on mass. FEBio's implicit solver stalls on exactly those.
- The source's AVW and PVW are two separate walls with no side connection; in FEBio we joined them into a tube with curved side
  bands ("curves", NOT IN SOURCE) and added a strap from the PM arch to the curve tip, contacts (curve self-contact, curve-strap,
  canal-PM arch) and some damping. Every such addition is marked "NOT IN SOURCE" in the run's `.feb.changes.txt`.
- The current blocker is the distal curve (next to the PeB): under load it folds into a hairpin and passes through itself and the
  strap (measured by `curve_cross.py` / `strap_cross.py`; under ~1 mm is acceptable), and late in the hold (t 1-1.5, load held
  constant) the canal could collapse.

## How we are fixing it (state on 2026-10-23)
1. **Canal width = source width** (L423b): our walls + curves had been 8-28 % wider than the source; narrowing them delayed the
   late collapse and changed the prolapse much.
2. **Softer strap** (strap material x1/3, L431): first run to finish t 1.5 cleanly with no late collapse (old wall material).
3. **Tissue that folds instead of resisting compression**: the source's Yeoh/Marlow wall law is stiffer in compression than in
   tension, so the wall gets pushed through itself instead of buckling. We are moving the curves and the AVW / PVW to a
   tension-only fibre material (soft matrix + in-plane fibres that only carry tension), fitted by a one-element strip test
   (`claude_diag/shell_tests/curve_fibre/gen.py`) to the source law in tension.
   - 2 fibres (0 / 90 deg) failed (L442c): no stiffness for a diagonal pull and too low a bulk modulus; the wall stretched and
     its elements changed volume until the run blew up.
   - 4 fibres (0 / 45 / 90 / 135 deg) + the source bulk modulus (L444 / L445): halves the volume drift but not to the source
     level yet; L447 tries a stiffer matrix (compression 6x softer than the source instead of 14x).
4. **Spring attachments**: the PM connectors (AVW -> top of the PM arch) now go on every node of the AVW's distal row, not on the
   curve corners (L445, L448).

5. **(10-24)** The middle law ("softer in compression, not too soft") runs on the strap model (L465, the new base) but did not
   stop the curve corner penetrating itself: that is the open problem now.

## Where to look
- `CLAUDE.md` (how a Claude Code session starts and the run rules), `HANDOFF.md` (whose turn, what is running, next).
- `claude_diag/README_2026-10-09.md` section 8: the run-by-run journal (newest at the end); `claude_diag/run_log_notes.csv`.
- `notes/MEMORY_NOTES.md`: the working rules (e.g. one change per variant, no shared cores, say "NOT IN SOURCE").
- Bases: `claude_diag/runs/<run>/<run>.feb` + `.feb.changes.txt` (only the ones worth building on; results files stay local).
