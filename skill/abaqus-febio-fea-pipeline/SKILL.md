---
name: abaqus-febio-fea-pipeline
description: Use when working with Abaqus (.inp) or FEBio (.feb) finite-element model files -- converting or comparing model versions, tracing a converted FEBio model's boundary conditions/contacts/materials back against an original Abaqus source, correcting a conversion to match the Abaqus source's intent, or diagnosing FEBio solver convergence failures (negative jacobians, NaN, non-convergence, "runs fine without part X but fails with it", "won't get past time/load Y"). Push to use this whenever the user mentions .feb or .inp files, FEBio, FEBioStudio, Abaqus finite-element/structural models, or a biomechanical/structural simulation that won't converge or run -- even if they don't say "FEA" or "finite element" explicitly. Also covers general techniques for safely editing and comparing multi-megabyte FE text/XML files.
---

# Abaqus <-> FEBio FEA pipeline

## Scope

This skill covers the *comparison, conversion-fidelity, and convergence
debugging* pipeline for Abaqus/FEBio structural FE models: diffing model
versions, cross-referencing a converted FEBio model against its original
Abaqus `.inp` source to find and fix conversion drift, and systematically
diagnosing why a FEBio model won't run to completion. One specific piece of
material-parameter fitting is covered as part of conversion fidelity:
re-fitting an Ogden hyperelastic law from an Abaqus `*Uniaxial Test Data`
table when the target solver doesn't support Marlow directly (`scripts/fit_ogden.pl`,
`abaqus-inp-format.md`, `febio-xml-format.md` gotcha 16). Broader inverse-FEA
/ parameter-optimization-against-simulation-output / PCA-ML workflows built
on top of this pipeline are **not** yet covered -- that is a related but
distinct extension, expected to be merged into a later version of this
skill from other reference conversations. If asked to do that kind of
broader work, treat the material below (including the Ogden-fitting piece)
as the foundation but don't assume the rest is covered.

Everything here generalizes across whatever specific anatomical/mechanical
parts, node ranges, or file names a given task involves -- the domain
knowledge that matters is about the **file formats and solver behavior**,
not any particular model's geometry.

## The workflow

Real tasks in this space are usually some subset of these five stages, in
roughly this order:

1. **Understand what changed.** Given two versions of a model (e.g. "what
   changed besides the thing I asked for?"), diff them section by section
   rather than as raw text -- huge XML/keyword files have enormous amounts
   of coordinate-precision noise that a line-diff drowns in. See
   `references/geometry-cross-referencing.md` Rule 1.

2. **Trace FEBio constructs back to the Abaqus source.** A converted model
   almost always has an original `.inp` sitting behind it, and it is the
   ground truth for what the model is *supposed* to represent -- element
   types, boundary conditions, tie/contact topology, material laws. Read
   `references/abaqus-inp-format.md` for the keyword structure, and match
   geometry between the two files **by coordinate, never by ID**
   (`references/geometry-cross-referencing.md` Rule 3) -- node numbering is
   never preserved across a conversion. Three source checks that are easy to
   skip and expensive to miss: the `*Step` procedure (an explicit-dynamic
   source tolerates near-mechanisms a static FEBio solve cannot), every
   `*Contact Pair` inside the step, and every `*Tie`/`CONN3D2` touching the
   part in question, resolved into actual node pairs. Converters routinely
   drop or partially convert these (`references/abaqus-inp-format.md`).

3. **Fix mismatches to match the Abaqus source's intent**, being explicit
   about which fixes are exact conversions (e.g. re-deriving the correct
   boundary-condition node sets) versus fidelity approximations that are
   unavoidable because the target format lacks a direct equivalent (e.g.
   refitting a Marlow hyperelastic law as Ogden when the target solver
   doesn't support Marlow -- `scripts/fit_ogden.pl`,
   `references/abaqus-inp-format.md`, `references/febio-xml-format.md`
   gotcha 16 -- or working around FEBio having no bare-node contact surface,
   `references/febio-xml-format.md` gotcha 5). Say out loud which is which;
   don't let an approximation pass as an exact match. If the underlying mesh itself is genuinely
   defective (a sliver/degenerate element baked into the geometry, not
   caused by deformation), see `references/mesh-repair-with-mmg.md` for
   boundary-preserving remeshing -- critically, preserve every original
   boundary node rather than letting a remesher regenerate the boundary
   too, or the "fixed" region can end up numerically stable but
   anatomically/structurally connected to the wrong place.

4. **Run it and read what FEBio itself says.** Command-line `febio4.exe`,
   not just the GUI -- see `references/convergence-debugging.md` for the
   log vocabulary and exact invocation. A model that doesn't run is not a
   dead end, it's a diagnostic starting point.

5. **If it doesn't converge, debug systematically**, one variable at a
   time, keeping a comparison table across attempts, rather than guessing.
   Full framework in `references/convergence-debugging.md`.

## Before touching a real multi-megabyte file, know these things

These are the highest-cost lessons from doing this work; full detail and
more gotchas are in the reference files, but these are worth holding
in mind constantly:

1. **FEBio requires node IDs to increase monotonically across the entire
   `<Mesh>` section**, not just within one `<Nodes>` block. Inserting new
   high-numbered nodes in the "logically right" spot in the middle of the
   file, ahead of pre-existing lower-numbered blocks, silently breaks
   parsing hundreds of lines later with a confusing error. New nodes
   belong at the end of all node declarations, after the current
   highest-ID block. -> `references/febio-xml-format.md` gotcha 1.

2. **Never guess unfamiliar FEBio XML syntax directly on the real file.**
   Write a 20-40 line throwaway `.feb` isolating just the construct in
   question and run it -- seconds per iteration, unambiguous errors,
   instead of a slow guess-and-check cycle against a huge file where the
   error's line number may not even point at the real cause. ->
   `references/febio-xml-format.md` gotcha 6 has a ready-to-use minimal
   skeleton.

3. **Match Abaqus and FEBio geometry by exact coordinate, never by node
   ID** -- Abaqus part-local IDs and FEBio's flat renumbered IDs have no
   relationship. -> `references/geometry-cross-referencing.md` Rule 3.

4. **A "more correct" fix that makes the model less stable is real
   information, not a contradiction.** If replacing an approximate but
   long-working configuration with a precise one causes an *earlier*
   failure, something else in the model was quietly depending on the
   imprecision -- isolate by testing old+new together before assuming the
   precise version is simply wrong. -> `references/lessons-learned.md`
   has this exact case.

5. **Check for the user's own prior attempts before re-deriving anything.**
   Differently-named file variants, stray `.log`s, anything that looks
   like earlier iteration -- these routinely contain hard-won parameter
   values that would otherwise cost another full debugging pass to
   rediscover. Run the actual current file fresh rather than trusting a
   log that might predate later edits. -> `references/geometry-cross-referencing.md`
   Rule 10.

6. **FEBio's `tied-node-on-facet` pulls nodes onto the surface; it doesn't
   keep the gap** the way Abaqus `adjust=no` does. For an Abaqus node-to-node
   tie, or a tie to truss/wire nodes, use `<Constraints>` linear constraints
   instead: exact, with no facets or phantom nodes (the phantom-ribbon
   workaround is singular by construction). -> `references/febio-xml-format.md`
   gotchas 17-18, 8.

7. **When the same few elements fail in every variant, check their undeformed
   geometry before changing any physics**, including through the thickness
   for shells. A thick shell with averaged nodal normals can be half-inverted
   before any load at mesh folds; `<shell_normal_nodal>0</shell_normal_nodal>`
   fixes it. Material, formulation and solver changes won't.
   -> `scripts/shell_pinch_check.py`, `febio-xml-format.md` gotcha 19.

8. **Mesh surgery (remeshing, re-export, node duplication) can silently orphan
   boundary conditions.** The elements move to new node IDs at the same
   coordinates, and the BC still pins the old ones. Audit active-versus-declared
   nodes per BC across versions after any such step; don't assume the
   "baseline that runs" is correctly supported. -> `scripts/bc_audit.py`,
   `geometry-cross-referencing.md` Rule 12.

9. **A converted `<NodeSet>` existing is not proof its BC survived conversion.**
   A converter can preserve a node grouping under the right name with the
   right members while dropping what used it -- no `<bc>` ever references the
   set. The model runs and converges normally for a while, so this looks
   exactly like a material/contact/mesh problem rather than a missing
   constraint. Cross-check every source `*Boundary` set against the target's
   actual BC usage (not just its node-set declarations) whenever a
   converged-then-stuck model doesn't have an obvious cause. ->
   `abaqus-inp-format.md`.

10. **A refit that matches the uniaxial curve can still be 2-3x too soft in
   biaxial stretch.** Abaqus Marlow is I1-only, so refit it to FEBio `Yeoh`
   (also I1-only) with `scripts/fit_yeoh.py`, not to a 1-term Ogden. Also
   choose `k` for the working strain range. A `k` from mu0 matches Poisson's
   ratio only at zero strain, and the material loses volume stiffness as it
   stiffens. Separately, a single shell element pulled only through its nodes
   reads 2-3x too soft because the back face lags. Prescribe
   `sx`/`sy` as well. -> `febio-xml-format.md` gotchas 23-24.

11. **A spring chain in compression buckles, and an explicit source rides through on mass.** When a wall
   starts at a node of a 1D spring/truss chain, check the sign of the chain's axial force. If the
   source's trusses are heavy (`*Density`) under `*Dynamic, Explicit`, give the FEBio chain that mass
   (mass-only truss elements, written as `line2` + `<BeamDomain type="linear-truss">`) and run DYNAMIC. A dynamic run that reaches t = 1 is still not an
   equilibrium: check nodal speeds at the end, and settle with held loads plus `mass damping`.
   -> `convergence-debugging.md` (diagnostic loop), `febio-xml-format.md` gotcha 24.
   Before that, check which way the chain gives way and whether a source support in that direction
   is missing. Item 12 has the case where it was.

12. **A loft can lose its anchor without looking disconnected, and it needs its own material.** Lofts
   (shell surfaces lofted over the source's connector lines and meshed with MMG; in the reference model
   the domains are named `*_fan`, but call them lofts) stand in for families of `CONN3D2` connectors.
   When the structure a loft was lofted onto is replaced (e.g. a tube by a truss chain), its edge can be
   left free while it still shares nodes with the tissue. Audit every connector family by both ends,
   and check what holds each loft's boundary nodes in the current version. **Whenever a loft is
   created or reattached, give it its own material fitted to the Abaqus elements it replaces** (strip
   model; a 1-term Ogden follows J-shaped connector curves). A shared generic material was 16-60x too
   stiff in the reference model. Also take k from the stiffened modulus (about 250 mu0, not 50 mu0) and use
   the other lofts' thickness, or the loft swells or wrinkles. Fitted this way, it matched the source's
   connectors like for like. **Then audit every loft, not only the one being debugged.** In the reference
   model the other 12 lofts all carried one material inherited from a rigid Abaqus beam section (E scaled
   down, but the beam's density, so the lofts outweighed all the tissue 24x). Against their own connectors
   they ranged from 0.6x (too soft) to 10^7x too stiff, the latter where the source connectors were
   impaired to almost nothing. Trust a strip fit to about +-50 % and check it in the running model: compare
   the connector-end elongations with a run that has the connectors as springs. That check found the fitted
   lofts right for most families, but no loft can follow connectors that shorten: a loft carries
   compression, and a connector table that starts at (0, 0) doesn't. The spring version also converged
   best, so it is a usable model in its own right.
   -> `abaqus-inp-format.md` (lofts and tubes, loft material fit, connector audit, in-situ check).

## Command-line basics

```bash
# Typical FEBioStudio install location on Windows; adjust if different
"/c/Program Files/FEBioStudio/bin/febio4.exe" -i model.feb
```

Run from a scratch directory per variant when comparing configurations --
FEBio writes `<name>.log` and `<name>.xplt` next to the input. If no `.log`
appears, the run failed to parse; read stdout for the `tag "..." (line N)`
message rather than looking for a log that was never written.

For per-element failure detail instead of just an aggregate count, run
`febio4.exe` with no arguments to get its interactive console, then
`set output_negative_jacobians 1` followed by `run model.feb` -- and pin it
single-threaded (`OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`) if the specific
failing element matters, since the default multi-threaded solver isn't
deterministic run-to-run. Full detail in `references/convergence-debugging.md`.
Non-interactive form for background runs:
`printf 'set output_negative_jacobians 1\nrun -i model.feb\nquit\n' | febio4.exe > model.console.txt`,
then `scripts/negjac_report.py` maps the reported elements to domains.

The Python scripts below need Python 3 + numpy. On Windows, `python` may be
the Microsoft Store stub; `py -0p` lists the real interpreters (e.g. an
Anaconda install used as `py -3.10`). Run them from a short path, since deeply
nested scratch folders exceed the 260-character path limit.

Validate XML well-formedness (cheap, catches a lot, but is not a substitute
for FEBio's own stricter semantic checks) before every run once you've hand-edited
the file:

```powershell
try { [xml]$xml = Get-Content -Raw "model.feb"; "XML VALID" }
catch { "XML ERROR: $($_.Exception.Message)" }
```

## Bundled scripts

Ready-to-use instead of re-deriving the same awk one-liners each time:

- **`scripts/check_monotonic_node_ids.sh model.feb`** -- flags every place a
  `.feb` file's node IDs decrease across `<Nodes>` blocks (gotcha 1 above).
  Run this on any file you've hand-edited, before running FEBio on it.
- **`scripts/build_node_lookup.sh model.feb > nodes.txt`** -- extracts a flat
  `id|x,y,z` table for fast coordinate cross-referencing.
- **`scripts/nearest_point_distances.sh targets.txt reference.txt`** -- for
  each point in one lookup table, the minimum distance to any point in
  another. Use it to measure a real contact gap before picking a tolerance,
  or to check how far a converted surface has drifted from an original
  Abaqus centerline.
- **`scripts/validate_feb_xml.ps1 -Path model.feb`** -- well-formedness check
  before every solver run on a hand-edited file.
- **`scripts/nearest_facet_distance.pl targets.txt nodes.txt facets.txt`** --
  like `nearest_point_distances.sh` but measures true point-to-*facet*
  (triangle) distance, not point-to-vertex -- use this specifically when
  calibrating a proximity-based tied contact's tolerance, since the real
  engagement gap can be meaningfully smaller than the nearest-vertex
  distance suggests (`geometry-cross-referencing.md` Rule 8).
- **`scripts/fit_ogden.pl data.txt [max_N]`** -- fits an N-term Ogden
  hyperelastic model (default tries N=1..3, reports all) to an Abaqus
  `*Uniaxial Test Data` table via multi-start coordinate-descent least
  squares, no Python/scipy required. Outputs the Abaqus/Ogden-1972
  `(mu_i, alpha_i)` and the converted FEBio `(c_i, m_i)` per term
  (`febio-xml-format.md` gotcha 16) plus R^2, a Drucker-stability flag, and
  a monotonicity flag (a stronger overfitting tell than R^2 alone) for each
  N -- pick the smallest N where R^2 has saturated, not the highest-R^2 one
  by default. For a *Marlow* source, prefer `fit_yeoh.py` (item 10).
- **`scripts/fit_yeoh.py data.txt [--nmax 4] [--nu 0.47] [--ogden c1,m1,...]`**
  -- Marlow `*Uniaxial Test Data` -> FEBio uncoupled `Yeoh` (I1-only, like
  Marlow). It is a linear least-squares fit that reports N=1..nmax, picks the
  smallest stable fit within 10 %, and prints the FEBio block. It also gives
  `k` from mu0 and the `k` that reproduces Abaqus' constant-nu volume change
  at 20-50 % strain. Then it compares uniaxial/planar/equibiaxial stress
  against the Marlow construction and an optional existing Ogden fit.
- **`scripts/feb_postmortem.py model.feb model.xplt [model.log]`** -- where a
  failed run deforms and blows up: per-domain displacement history,
  relative-volume (J) extremes with element IDs, and, if the failed
  iterations were plotted, which nodes and domains the divergence starts in.
- **`scripts/negjac_report.py model.feb console.txt`** -- turns the console's
  per-element negative-jacobian reports into a per-domain/per-element table.
- **`scripts/shell_pinch_check.py model.feb --domains A,B`** -- flags shell
  elements whose offset faces are pinched or inverted at rest along averaged
  nodal normals (the thick-shell defect in item 7).
- **`scripts/bc_audit.py v1.feb v2.feb ...`** -- active-versus-declared nodes for
  every BC across versions, plus the live nodes an orphaned BC should
  probably pin (item 8).
- **`scripts/summarize_feb.py model.feb [--collapse PREFIX]`** -- one-line-per-item
  structural summary of a `.feb` for diffing versions section by section
  (Rule 1 of the geometry reference), with mesh blocks reduced to counts and hashes.
- **`scripts/xplt_reader.py`, `scripts/feb_model.py`** -- the plot-file reader
  (safe on multi-GB and still-being-written files) and the `.feb` parser the
  other scripts import. They are usable directly for custom checks.

## Reference files

Read these as needed -- they're where the deep, hard-won detail lives:

- **`references/abaqus-inp-format.md`** -- Abaqus keyword structure, what
  matters for extraction/comparison, extraction one-liners, the specific
  quirks that bite during a conversion (part-local node IDs, `generate`
  range semantics, `*Tie` surface types, `PINNED`/`ENCASTRE`/`XSYMM`
  boundary condition semantics, Marlow-hyperelastic-on-a-truss reducing to
  a plain force-strain curve), plus resolving `*Tie` slave/master pairs,
  `CONN3D2` connectors (force-first tables, extrapolation, counting them
  per part), assembly-level reference nodes and instance-less `*Nset`s, the
  `*Step` checks (explicit dynamics, `*Contact Pair`), and a `*Boundary` set
  that survives conversion as a correctly-populated FEBio `<NodeSet>` with no
  `<bc>` ever referencing it -- a silent, graceful-looking conversion gap,
  distinct from the orphaned-BC case in item 8 below. Also how the reference
  project's lofts and tubes were made (lofts over connector lines, meshed
  with MMG; tubes as separate connection surfaces), what a loft needs
  (anchored at both ends, its own material fitted to the connectors it
  replaces), and auditing connector families by both ends (item 12).
- **`references/febio-xml-format.md`** -- FEBio's section layout and
  twenty-five specific XML/semantic gotchas, roughly ordered by how much time
  each one cost, including the node-ID-monotonicity rule, where `DiscreteSet`
  vs. `discrete_material` actually live, why a `<Surface>` can't be built from
  bare nodes, how to extract valid type-name keywords straight out of
  the installed `.dll`s (and, when that's not enough, getting a real
  working `.feb` exported from FEBioStudio as unambiguous ground truth for
  exact syntax; also why a same-sounding element formulation like `udg-hex`
  can still be far stiffer than the Abaqus element it looks like it should
  replace), the three tied-contact types and how they differ
  (including which one has no `auto_penalty` and needs its `penalty` scaled
  to your own model's materials, not copied from an example), the
  `max_distance="0"` = unlimited trap, the multi-target "double-tie" bug,
  facet-winding-consistency propagation, a confirmed-working discrete-spring
  tie pattern (distinct from a similar-looking but different failure mode),
  why fixing a shell's nodes doesn't fully immobilize it, the
  Abaqus<->FEBio Ogden parameter conversion, `tied-node-on-facet` closing
  initial gaps, exact ties via `<Constraints>` linear constraints, thick-shell
  pinching and `shell_normal_nodal` (including the correct FEBio-4 top-face
  offset geometry, not the older mid-surface one), the available shell
  formulations, per-iteration plot output, the `<offset>` a shell contact
  needs when the Abaqus source was meshed mid-surface, and why
  `<discrete dmat="k">` resolves by list position, not by matching `id`.
  Gotcha 23 covers refitting a Marlow material as I1-only `Yeoh` rather than
  uniaxial Ogden, and choosing `k` beyond mu0. Gotcha 24 covers
  single-element shell tests needing back-face (`sx`/`sy`) BCs, the verified
  shell-displacement BC syntax, and a delayed smooth-step load curve.
  Gotcha 25 covers FEBio 4.13 beams (`line2` + `<BeamDomain type="linear-beam">`,
  shear locking, the two-pin mechanism). A beam used as a stabiliser on a tied
  spring chain made the reference model fail earlier in every variant tried.
- **`references/geometry-cross-referencing.md`** -- techniques for working
  with huge FE text files without loading them whole: ID->coordinate lookup
  tables, coordinate-based cross-referencing between Abaqus and FEBio,
  distinguishing a genuinely merged/welded node from a merely coincident
  one, mesh-quality checks without a GUI, safe large-block editing (and
  the reverse-line-order trap), the facet-only "phantom surface"
  pattern for converting a bare 1D chain into something FEBio's
  contact/tie machinery can use, true mesh-boundary-edge detection (vs.
  the much easier-to-mistakenly-use "not shared with anything" check),
  correctly extracting a solid element's boundary faces (the hex8
  node-index table facet-winding bugs usually trace back to), and
  measuring a tied-contact's real point-to-facet gap rather than
  point-to-vertex, checking a thick shell's through-thickness geometry
  (Rule 11), and auditing boundary conditions after mesh surgery (Rule 12).
  Rule 6 (phantom surface) is kept for history but superseded.
- **`references/convergence-debugging.md`** -- the full systematic
  diagnostic loop: how to read a FEBio log, getting per-element
  negative-jacobian detail via the interactive console (and why that
  diagnostic run must be pinned single-threaded -- FEBio's default
  multi-threaded solver is not deterministic run-to-run), the diagnostic
  ladder of things to try (isolate new contacts -> trace indirect shared-
  geometry load paths -> check engagement -> ramp via a load curve -> expect
  contact-parameter sweeps to be chaotic/non-monotonic -> vary solver
  patience, including why `max_retries`x`max_refs` multiply and when
  `lstol=0` helps -> full Newton -> arc-length -> dynamic analysis ->
  stiffness scaling as a buckling-vs-mechanism diagnostic -> check actual
  load fraction reached), and when to stop iterating blindly and ask for
  human eyes on the deformed shape instead. Also: localizing a failure from
  the plot file's failed iterations, the decouple test for "runs without
  part X", distinguishing global softening from a local event, checking the
  source's analysis type first, and running single-variable variants in
  parallel.
- **`references/lessons-learned.md`** -- a quick-scan table of specific
  mistakes made, how each was caught, and the generalizable lesson from
  each -- read this first if you want the fastest orientation on what
  *not* to repeat.
- **`references/mesh-repair-with-mmg.md`** -- boundary-preserving surface
  remeshing with `mmgs` for a genuinely degenerate/sliver mesh region:
  getting `mmgs` without admin rights, the `.mesh` format, deriving
  `hmax`/`hmin`/`hausd` from the patch's own characteristic edge length
  instead of guessing, and why regenerating the boundary along with the
  interior can produce a numerically-stable but anatomically-wrong result.

## A note on rigor

Nearly every fix in the reference session that "should" have worked based
on general FE intuition was tested empirically rather than trusted on
reasoning alone -- and several confident-sounding hypotheses turned out to
be backwards (a finer solver step size performing *worse*; a more
mechanically-faithful reconstruction being *less* stable than the
approximation it replaced). Treat every proposed fix in this domain as a
hypothesis to A/B test against a recorded baseline, not a conclusion to
apply and move on from.
