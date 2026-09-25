# Systematically debugging FEBio convergence failures

The goal of this process is to turn "it doesn't converge" into a specific,
falsifiable diagnosis before trying to fix anything. Changing parameters
without first characterizing the failure produces a lot of wasted solver runs
and conclusions that don't generalize. Everything here assumes command-line
`febio4.exe`, not the FEBioStudio GUI -- the GUI is fine for building/visually
inspecting a model, but the command line is far faster for the run-edit-run
loop this process needs.

## Running FEBio from the command line

```bash
"/c/Program Files/FEBioStudio/bin/febio4.exe" -i model.feb
```

Run it from (or `cd` into) the directory containing the model, or copy the
`.feb` into a scratch run directory first -- FEBio writes `<name>.log` and
`<name>.xplt` next to the input file, and you'll want a clean directory per
variant when comparing many configurations side by side. For a long-running
model, launch it as a background task and poll the `.log` file rather than
blocking on the tool call.

**If no `.log` file appears at all**, the run failed before logging even
started -- almost always an XML parse error. Check the raw stdout/stderr for
`Reading file ... FAILED!` and the specific `tag "..." (line N) : ...`
message (see `febio-xml-format.md` for the common causes). Don't spend time
looking inside a log file that was never created.

A different, unrelated cause of instant failure with no `.log` at all: running
`febio4.exe` from a working directory with a very long path (deeply-nested
temp/scratch folders are the usual culprit) can hit a Windows path-length
limit and fail before even reading the input file, with an error like `path
longer than allowed for a Win32 working directory` or the shell simply
reporting the command as not found. This has nothing to do with the model
file's own content -- `cd` somewhere short before running, or copy the
`.feb` to a short-path scratch directory first.

## Reading the log: the vocabulary

```
Number of time steps completed .................... : N
```
The single most important line for "how far did it get" -- compare this (and
the corresponding `time=` value on the last successful step) across variants
as the primary success metric, not the number of solver iterations.

```
Nonlinear solution status: time= T
   convergence norms :     INITIAL         CURRENT         REQUIRED
      residual            ...
      energy              ...
      displacement        ...
```
`INITIAL` is the residual right after the load/BC increment is applied but
before any Newton correction; `CURRENT` is after the latest iteration;
`REQUIRED` is the convergence tolerance (a fraction of `INITIAL`, controlled
by the solver's `dtol`/`etol`/`rtol` settings). A huge, unphysical `INITIAL`
displacement norm (tens to hundreds of length units, when the model itself
is a few tens of units across) at the very first iteration of a newly-active
contact/tie is a strong signal that a real geometric gap is being closed
instantly ("snapping") rather than gradually -- see the ramping technique
below.

```
WARNING: Problem is diverging. Stiffness matrix will now be reformed
ERROR: N negative jacobians detected.
ERROR: NAN detected
```
These are different severities and point in different directions:
- A **negative jacobian** after normal-looking residual convergence
  attempts is consistent with genuine large deformation/element inversion --
  investigate whether the mesh is being pushed through a real geometric
  limit (buckling/snap-through) or a self-intersection.
- **NaN**, especially appearing at a vanishingly small load fraction (well
  under 1% of the intended final load) or even with a near-zero contact
  penalty already ramped down, is much more consistent with a genuine setup
  bug -- a singular/ill-conditioned DOF, a degenerate (zero-length or
  zero-area) element, or a unit/scale mismatch -- rather than a real physical
  instability. Don't reach for "needs more solver patience" fixes (smaller
  steps, more retries) for a NaN; they rarely help and can make interim
  diagnosis slower.

```
WARNING: No contact pairs found for tied interface "X"
WARNING: N isolated vertices removed
```
See `febio-xml-format.md` gotchas 9-10 -- the first means a contact is
silently doing nothing (check the real gap vs. the tolerance), the second is
usually harmless orphaned geometry.

## Getting per-element detail: the interactive console, and why it must be single-threaded

The plain `-i` invocation only ever reports an aggregate count ("N negative
jacobians detected") -- it never says *which* element. To get the specific
element ID (and gauss point), run `febio4.exe` with no arguments at all to
drop into its interactive console, then issue two commands before running:

```
febio4.exe
febio> set output_negative_jacobians 1
febio> run model.feb
```

This makes every negative-jacobian report print as `Negative jacobian was
detected at element N at gauss point G` instead of just a count, which is
usually the fastest way to go from "something is failing" to "this specific
element/domain is failing" -- cross-reference the element ID against each
domain's `<Elements>` block to find which named domain it belongs to.

**Critical caveat: FEBio's default multi-threaded (pardiso) linear solver is
not deterministic run-to-run.** Thread-scheduling-dependent floating-point
summation order means two runs of the *identical* file can converge via
different iteration paths and fail at a different time step -- sometimes by
a large margin. This means a diagnostic run (via the interactive console
above) can silently diagnose a *different* failure than the one actually
being investigated, unless both are pinned to a single thread:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 febio4.exe -i model.feb
```

Verify determinism is actually in effect for a given install by running the
same file twice with these variables set and diffing the two logs -- they
should be byte-identical. Once confirmed, **always set these two variables
before any run whose specific failure point or failing element matters** --
which is effectively every diagnostic run. For a run where only "does this
get further, yes or no" matters (a quick A/B check, not a diagnosis),
multi-threading is fine and considerably faster (roughly 4-8x) -- reserve
single-threading for the run you're actually going to read element IDs out
of. Single-threading also makes each stuck time step much more expensive in
wall-clock terms (no parallel speedup on top of however many retries/
reformations it needs), so a diagnostic run on a large model can take
substantially longer than the equivalent multi-threaded check -- factor that
into how long to wait before concluding a run is stuck rather than just slow.

Non-interactive form (works as a background job; both `run model.feb` and
`run -i model.feb` are accepted), plus a per-domain summary of the reports:

```bash
printf 'set output_negative_jacobians 1\nrun -i model.feb\nquit\n' | \
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 febio4.exe > model.console.txt
py -3 scripts/negjac_report.py model.feb model.console.txt
```

## Localizing a failure from the plot file

The console tells you which elements went negative; the plot file tells you
where the model *started* going wrong. Plot the failed iterations too
(`<plot_level>PLOT_MINOR_ITRS</plot_level>`; FEBioStudio's debug runs already do
this, and it can take GBs, see `febio-xml-format.md` gotcha 20), then:

```bash
py -3 scripts/feb_postmortem.py model.feb model.xplt [model.log]
```

It prints each domain's displacement history, the J (relative volume)
extremes per domain at the last converged state with element IDs, and, for
each failed iteration, which nodes' increments blow up and which domains they
belong to. In the reference session the log said only "2 negative jacobians
detected". The post-mortem showed every failed iteration diverging at the
junction of one muscle shell and two fans, which is what pointed at the
pinched thick shell (`febio-xml-format.md` gotcha 19).

Two readings of that output are worth distinguishing:
- **Global softening vs a local event.** Compare how compliant the model is
  across variants at the same load: (first-iteration displacement norm)^0.5 /
  (dt x load-curve slope), per converged step. If a variant's compliance
  tracks the baseline almost exactly and it still fails much earlier (it did in
  the reference session: within 1-3 % right up to the failure), the failure is
  a local event, and globally stiffening or softening anything will not help.
- **"Small strains everywhere" at failure.** If every domain's J is within a
  few % of 1 and displacements are a small fraction of the model size when
  it fails, suspect a setup/geometry defect (a near-mechanism, a pinched
  element, an unsupported DOF), not a material or large-deformation limit.

## The diagnostic loop

0. **Check what analysis the source model ran** (`abaqus-inp-format.md`,
   "Check the `*Step`"). An Abaqus `*Dynamic, Explicit` source tolerates
   floppy membranes, truss chains with free interior nodes and dangling flaps
   that a static implicit FEBio solve cannot. Knowing that up front tells you
   to hunt for near-mechanisms and degenerate reference geometry rather than
   to tune the solver.

1. **Get a baseline number.** Run the model as given (or the last known
   version) and record: time steps completed, the `time=` value reached, and
   whether the failure was immediate (0 steps) or gradual (several steps
   before failing). This distinction matters enormously: immediate failure
   at essentially zero load is a setup/configuration problem; gradual
   failure after meaningful load has been applied is more likely a genuine
   structural/numerical limit that needs a different technique, not just a
   bug fix.

2. **Change exactly one thing per run**, and keep every variant as a
   separate, clearly-named file/run-directory rather than overwriting --
   building a small comparison table (variant -> steps completed -> time
   reached -> nature of failure) across many single-variable changes is what
   actually produces a diagnosis, not any single change in isolation.
   Candidates, in roughly the order they're worth trying:
   - **The decouple test, for "runs without part X but not with it".**
     Keep part X in the model but disconnect it: give it its own copies of
     the nodes it shares with everything else (and drop its ties/connectors),
     so it is present but carries no load. If the run then matches the
     without-X baseline exactly (it did in the reference session: 0.6043 vs
     0.6043), X's body, material and DOF count are fine. The problem is how X
     is *coupled and loaded*, so test coupling variants next, not X itself.
   - **If the same few elements fail in every variant, inspect their
     reference geometry before changing any physics.** Include the
     through-thickness geometry for shells (`geometry-cross-referencing.md`
     Rule 11) and the BC audit after any mesh surgery (Rule 12). In the
     reference session, four physics-level changes (bulk modulus, shell
     formulation, thinner equivalent shell, dynamics) all stalled at the same
     place. The cause was elements that were already half-inverted before any
     load, and `shell_normal_nodal=0` fixed it.
   - **Isolate new/changed contacts.** If a change added or modified a
     tied/sliding contact, temporarily neutralize it (set its proximity
     tolerance to a value far smaller than any real gap, e.g. `1e-5`, so
     FEBio reports "no contact pairs found" and it becomes a no-op) and
     re-run. If that alone fixes or meaningfully improves things, the
     contact -- or specifically its *engagement*, not necessarily its mere
     presence -- is implicated. If several new contacts exist, disable them
     one at a time (or in combinations) to find which one(s) matter. Two
     ways to neutralize, useful for different questions: shrinking the
     proximity tolerance to a no-op value (above) tests whether the
     contact's *engagement* matters while leaving its `<SurfacePair>`/
     equations in the model; fully deleting the `<contact>` block plus its
     `<Surface>`/`<SurfacePair>` entries tests whether the contact's mere
     *presence* (extra equations/DOFs) matters at all, independent of
     whether it engages anything -- run both if the first test is
     ambiguous, since they isolate different variables.
   - **Trace indirect load paths through shared geometry**, not just
     direct connections, when a *newly added* contact/tie destabilizes an
     *unrelated, previously-stable* contact elsewhere in the model with no
     obvious direct link. Two contacts/ties can share load through a chain
     of ordinary shared-node connections even when neither directly
     touches the other: e.g. new-contact -> shares nodes with -> domain B ->
     shares nodes with -> the same solid that an already-validated,
     unrelated contact ties to. Loading that solid through the new path
     can push the old contact past its margin even though nothing about
     the old contact itself changed. Find this by checking, for each
     domain the new addition touches, which *other* domains it shares
     nodes with (not just its own direct contact partners), and following
     that chain until it reaches something the failing contact also
     touches. This is a real, sufficient explanation on its own -- but
     confirm it's not *also* something else by testing whether adjusting
     the newly-loaded contact's own penalty in isolation actually resolves
     the downstream failure; if it doesn't, the shared-load-path
     explanation may only be a contributing factor, not the whole story.
   - **Check whether a contact is even engaging.** Grep for "No contact
     pairs found". If a tie is supposed to be active but isn't, that's a
     tolerance/gap mismatch, not a stability problem -- measure the real
     gap (Rule 2/3 in `geometry-cross-referencing.md`) and set the
     tolerance to comfortably exceed it without being so loose it could
     match the wrong facet on a curved/folded surface.
   - **Ramp a newly-engaging contact in via a load curve**, rather than
     letting it apply full penalty force from the first increment. FEBio
     contact parameters like `penalty` accept an `lc="<id>"` attribute
     just like a pressure load does -- reference a load curve that ramps
     from 0 to 1 over the first few percent of the analysis:
     ```xml
     <penalty lc="2">0.0005</penalty>
     ```
     ```xml
     <load_controller id="2" name="TieRamp" type="loadcurve">
       <interpolate>SMOOTH STEP</interpolate>
       <extend>CONSTANT</extend>
       <points><pt>0,0</pt><pt>0.1,1</pt></points>
     </load_controller>
     ```
     This reduces (but does not always eliminate -- see below) the violence
     of an initial gap-closing snap. If the failure recurs later, once the
     ramp reaches full strength, that shows the *steady-state* engagement
     is the harder problem, not just the instant of first contact.
   - **Expect contact `penalty`/`max_distance` sweeps to be chaotic, not
     smooth.** Unlike most solver-patience settings, a proximity-based
     contact's `penalty` or `max_distance` can be *non-monotonic*: a value
     can outperform both a lower and a higher neighbor (e.g. `0.0005`
     converging further than either `0.0003` or `0.0007`), and a threshold
     that separates "works" from "fails" can do so while engaging the
     *exact same set* of nodes/facets on both sides -- meaning the
     parameter's raw numeric value is itself feeding some internal
     proximity-weighting/kernel near the cutoff, not just acting as a hard
     include/exclude gate. Confirm this distinction (same engaged set vs.
     different) with the true point-to-facet distance, not point-to-vertex
     (`scripts/nearest_facet_distance.pl` -- a facet's closest point can be
     meaningfully nearer than any of its own vertices, so a naive
     nearest-vertex measurement can overstate the real gap by a wide
     margin). If a binary search finds a boundary where the engaged set is
     provably identical on both sides, that rules out "a specific bad
     point crossed the threshold" as the explanation and means further
     narrowing of that exact parameter is chasing numerical noise, not a
     real defect -- settle on a validated working value inside the stable
     region (not right at the edge of it) rather than continuing to
     micro-tune, and treat a sweep's local optimum as directly tied to the
     current model/mesh, not something safe to reuse elsewhere without
     re-validating.
   - **Vary solver patience**, but treat this as a secondary lever, not a
     first resort: `max_retries`, `dtmin`, `max_refs`, `cutback` in the
     time-stepper/solver block. **Check that `cutback` is live:** FEBio reads it
     only with `<aggressiveness>1</aggressiveness>`. With the default 0, each retry
     takes dt0/(max_retries + 1) off the step, a linear walk down (FECore
     `FETimeStepController::Retry()`), so a `cutback` in the file may never have
     run. The log shows which rule is active: failed-attempt times that step down by
     a constant (0.69926, 0.699057, 0.698854 ... = dt0/21 with max_retries 20) are the
     linear rule. In the reference project every file had carried an inert
     `cutback` from the start. Counterintuitively, more patience
     (smaller minimum step, more retries, more stiffness reformations)
     sometimes makes a specific failure *worse*, not better, especially
     once other structural fixes have already changed the character of the
     problem -- always re-test whether a "more patient" configuration still
     helps after you've changed something else, don't assume it composes.
     `max_retries` (time-stepper: how many times a stuck step gets
     re-attempted at a smaller dt) and `max_refs` (solver: how many
     stiffness reformations one single attempt is allowed before it's
     itself declared failed) **multiply**, not add -- a stuck step's
     worst-case cost is roughly `max_retries x max_refs` full nonlinear
     iterations, each a real linear solve. On a large model this can turn
     one stuck step into a genuinely multi-hour wait. Lowering both
     together (e.g. `max_retries` 20->8, `max_refs` 25->12) is a legitimate
     way to get faster iteration during debugging even when it isn't
     expected to change *whether* something ultimately converges -- just
     re-verify against the original values once a configuration looks
     good, in case the tighter caps were quietly cutting off an otherwise-
     recoverable retry. If the specific symptom is a repeated `WARNING:
     Zero linestep size. Stiffness matrix will now be reformed` followed
     eventually by `Max nr of reformations reached` (the line search
     keeps finding no improving direction, not a geometry/jacobian
     problem), try `<lstol>0</lstol>` in the solver block to skip the
     line-search convergence requirement and accept the Newton step
     directly -- this trades away some protection against overshoot on
     large steps elsewhere in the model, so treat it as targeted for that
     specific symptom rather than a default to leave on everywhere.
     **Re-test it when the model changes a lot.** In the reference project
     `lstol` 0 cured such a stall in an early model version and then stayed
     in every later file. A small folding test (three thin shells, two of
     them folding, tied node-on-facet; full Newton, `lstol` 0.9) ran
     cleanly. One change each: line search off -> failed at 46 % of load.
     The project's whole solver block (BFGS `max_ups` 10, `lstol` 0,
     `opt_iter` 25, `max_refs` 25) -> failed at 4 %. BFGS with the line
     search on still finished (11 failed attempts vs 0). For folding shells
     the line search is the ingredient that matters. It carried over to the
     full model (implicit dynamic, 54k equations). With the line search on as
     the only change, it reached full load with 6 failed attempts in 23 min,
     against 99 in 54 min, and gave the same solution to 0.1 mm. That also
     rescued a variant that had failed at 27 % of load. Full Newton alone, and
     the test's whole block, were slower. With `max_refs` 100 each failed
     full-Newton attempt took minutes, and the run stopped advancing in
     wall-clock terms.
   - **Try full Newton instead of quasi-Newton** (`max_ups=0` under
     `qn_method` forces a full stiffness reformation every iteration
     instead of BFGS updates) if the failure looks like Newton's method
     losing the plot near a difficult point. In practice this often
     reproduces the exact same failure point as quasi-Newton, which is
     itself useful evidence: it means the problem isn't really about the
     iterative update scheme.
   - **Try arc-length continuation** (`<arc_length>1</arc_length>`,
     `<arc_length_scale>`) if the failure pattern looks like a genuine
     snap-through/limit point under load control (large iterative
     corrections needed right at a specific load level, under a pure
     load -- not displacement -- boundary condition). This is the
     textbook-correct tool for that specific failure mode, but tuning
     `arc_length_scale` without documentation is largely guesswork; if an
     untuned guess (e.g. leaving it at 0 for "auto", or trying a value of
     1) doesn't work within one or two tries, it's often more
     time-effective to fall back to other techniques than to keep
     guessing blindly.
   - **Try dynamic instead of static analysis** as a way to add inertial
     stabilization through a difficult transient -- but verify first that
     every material actually assigned to an element has a nonzero density
     (a dynamic analysis needs real mass), and be aware this can also make
     things worse if the transient dynamics themselves are numerically
     stiff. Reference result: even with the source model's own densities
     (the Abaqus original was explicit dynamic), FEBio implicit `DYNAMIC`
     with generalized-alpha `rhoi=0.5` stalled at t~0.07, far *earlier*
     than the static runs, both with and without the problem part. Don't
     expect dynamics to rescue a static model that has near-mechanisms or
     pinched elements; remove those instead.
   - **Scale material stiffness up and down** (e.g. x0.1 and x0.01) purely
     as a *diagnostic*, not a proposed fix, to distinguish elastic
     buckling from a kinematic mechanism. If a large stiffness change
     (orders of magnitude) only shifts the failure point by a small
     fraction, the structure is behaving like an under-constrained
     mechanism (a near-zero-energy motion path that barely depends on
     material stiffness) rather than undergoing real elastic buckling
     (which would be much more sensitive to stiffness). This distinction
     tells you whether to look for a missing/misplaced constraint
     (mechanism) or accept a genuine physical limit and reach for
     continuation methods (buckling).
   - **Isolate one suspect element into its own domain and material,
     stiffened, as a targeted diagnostic** when a specific oddly-shaped
     element (a mesh sliver, a collapsed edge) recurs at every failure but
     you don't want to change anything else about the model. Move just that
     element into a new `<Elements>`/`<SolidDomain>` (or `<ShellDomain>`)
     bound to a duplicated material with every stiffness parameter
     multiplied by a large factor (e.g. 10x), and re-run. If the wall
     doesn't move, that element is not the limiting factor, even though it
     looked like the worst-shaped one in the model. Confirmed case: a hex
     with a 0.1 mm collapsed edge against 0.9-1.9 mm neighbors (scaled
     corner Jacobian 0.008, the worst in its domain) sat right on a contact
     surface that kept failing at the same time step across many unrelated
     variants -- 10x stiffer left the wall exactly where it was, correctly
     ruling the element out and redirecting the search to a genuinely
     untested part of the model (which turned out to be the actual cause).
   - **Check what fraction of the intended load was actually reached.**
     A load ramped by a smooth-step curve from 0 to 1 does not apply load
     linearly with the analysis's internal `time` value -- compute the
     actual curve value at the failure time (e.g. a cubic smoothstep
     `3t^2-2t^3`) before concluding "it got to 20% of the run" means "it
     survived 20% of the load." Failing at a tiny fraction of true load is
     a much stronger signal of a setup/support problem than failing near
     full load.

3. **When a change helps, try to combine it with other independently
   helpful changes** -- but verify the combination actually stacks (it
   often doesn't; two independently-helpful changes can be neutral or even
   worse together, especially solver-patience settings combined with a
   change that already fixed the underlying issue). Don't assume additivity.

   **Run the single-variable variants in parallel.** Build each variant
   programmatically from one base file (Python `xml.etree.ElementTree` works
   on a multi-MB `.feb`: new nodes go in a new `<Nodes>` block right after the
   max-ID block, DiscreteSets at the end of `<Mesh>`, `<Constraints>` after
   `<Contact>`). Write a `changes.txt` next to each one, give each its own
   directory, and launch them together with a few threads each
   (`OMP_NUM_THREADS=4 MKL_NUM_THREADS=4`, e.g. 5-8 runs on a 32-thread
   machine). A progress table of `grep "converged at time"` per log turns an
   afternoon of serial guessing into one batch. Multi-threaded runs are not
   bit-reproducible, but in the reference session repeated failure times agreed
   to about 0.3 %, which is fine for "how far does it get" comparisons. Pin
   one thread only for the run you read element IDs out of. Two independent
   defects can each look like "no effect" when fixed alone (a shell fix alone:
   0.20; a coupling fix alone: 0.20-0.27; both together: 0.56), so when
   isolated fixes disappoint, also test them combined.

   **When the failure starts at a node of a 1D spring/truss chain, check the
   sign of the chain's axial force before calling it "slack" or "a
   mechanism".** A chain of axial springs has no bending stiffness. In
   tension its sideways stiffness is the geometric term F/L. In compression
   that term is negative, so the chain buckles sideways and a static solve
   stalls at the buckling point. The residual is small, but the displacement
   corrections wander at the chain node and the nodes tied to it. In the
   reference model's Load-LA wall, both posterior-arcus chains were in
   compression over most of their length (-1 % to -5.7 % strain, up to
   -0.8 N). The most compressed springs were the two either side of the node
   where every failed attempt diverged, and interior nodes had already drifted
   1-2.5 mm off the chain line. Compute each spring's strain from the plot
   file at the last converged states: `claude_diag`-style `chain_force.py`
   using `xplt_reader.py`. Check the source too. The Abaqus trusses there
   carried a density about 10^5 times the tissue's (`*Density 0.00011`). Under
   `*Dynamic, Explicit` that mass carries the chain through the instability,
   but the FEBio conversion's springs are massless. So an implicit-dynamic
   retry only mimics the source if the chain nodes get that mass as well.
   In the reference model that is what got through. The results were:
   - implicit DYNAMIC with massless chains crawled at t = 0.25;
   - with the Abaqus truss mass (mass-only truss elements, `febio-xml-format.md` gotcha 24) it passed
     the static wall of 0.29-0.33 and ran to t = 1.0, but slowly (1384 retries, ~1.5 h on 4 threads);
   - with that mass x10 it ran to t = 1.0 in 58 retries (~30 min).
   A 10x stiffer (non-faithful) material on the buckling-driving part only got the static run from 0.29
   to 0.65.

   **A dynamic run that reaches t = 1 is not automatically an equilibrium.** Check nodal speeds
   (du/dt between the last converged states) at the end, and compare the dynamic displacement against
   a static run at a time both reach. Heavy parts make the dynamic path lag the static equilibrium. In the
   reference model the LA, which is tied to the heavy chains, was at 9 mm against 22 mm static at t = 0.29
   even with the source's own mass. At t = 1 it was still moving at ~180 mm/s. An explicit source with
   the same heavy parts lags the same way, so "the FEBio result bulges more than the Abaqus result" can
   be the Abaqus result not having settled. To get the equilibrium, run on past the ramp with loads held
   and mass damping switched on (dynamic relaxation). The masses change the path, not the
   equilibrium. In the reference model this settled the LA at about 59 mm within 0.1 s of hold, with speeds
   falling 180 -> 67 -> 43 mm/s. The compressed chains had flipped into tension (+31 to +56 %): the model
   had passed through the buckling into a stable hanging state. The settle phase then crawled in tiny
   steps (C = 20/s, 1.12 -> 1.15 in 1.5 h), so read the equilibrium from the first ~0.1 s of hold rather
   than waiting for t_end. Better, extrapolate it: fit y = a + b exp(-(t - t0)/tau) to the
   displacement percentiles (median / p90 / max) over the hold. A later settle crawled at 0.01 of time
   per 1.5 h near t = 1.66. The fit (rms 0.02-0.05 mm, tau 0.2-0.33 s) put the equilibrium 0.5-1 mm below
   the last state, and running on to t = 2 would have moved it only ~0.5 mm. So stop such a run, and
   report the fitted asymptote as an estimate.

   **Mid-ramp states lag too, and a slower ramp shows by how much.** The same model with its load ramped
   over 3 s instead of 1 s (step size scaled with it) had an LA displacement 2.5x larger at the same load
   level (65 %), then snapped through at 81 %. So a 1 s dynamic run's mid-ramp states are not points on the
   load-displacement curve. A slower ramp cuts speeds with the ramp time and inertial forces with its
   square: the cleaner route to the loaded state, since it adds no dissipation. It is not like for like with a
   source's own ramp timing, so keep the source's ramp for a time-history comparison. Mass damping from
   t = 0 (C = 20/s) cut the early failures but stalled at the same point. It also slowed the transient
   visibly, because the wall moved at 100-200 mm/s, ~10x the quasi-static speed the damping had been
   sized for. Measure the speeds before sizing damping.

   **Screening runs stalled where the damping switched on.** Six variants of a later model reached full
   load (t = 1). All then stalled within t = 1.00-1.05, where the mass damping ramped from 0 to C = 20/s in
   0.05 s. For screening, end at t = 1 and keep the settle for the chosen model. Otherwise ramp the damping
   on more gently (untested). Also run the control in the same batch: in two batches the control stalled
   (0.88, 0.92) where its own earlier copy had passed t = 1 easily, so this base is path-dependent. Read
   single-run differences as noise unless they are large. Test other C values or a static restart from the settled state as separate
   variants. The unfaithful alternative did not do as well: with the same dynamics, the source's material
   refit as a 1-term Ogden failed at 0.69, while the faithful I1-only law went to full load.

   **Before adding stiffness or mass, look for the sideways support the chain has lost.** Take each chain
   node's zig-zag drift (its displacement minus the mean of its two neighbours', across the chain). Split
   it along the candidate support directions: the normal of the sheet tied to the chain, and the direction
   of any source connectors or fans that should hold that node (`claude_diag`-style `chain_buckle_dir.py`).
   In the reference model the drift was entirely out of the LA's plane (1.01 of 1.01 mm at the failing
   node). The source's 26 arcus-to-vaginal-wall connectors pointed within 22-58 deg of it, and they were
   missing from the FEBio model (`abaqus-inp-format.md`, fans that lost their anchor). What a support has
   to supply: a compressed chain removes about 2F/L (one node) to 4F/L (alternating zig-zag) of sideways
   stiffness per node, 0.4-0.7 N/mm at F = 0.8 N and L = 4.2 mm. Results, one change each on the static
   run that stalled at t = 0.29 (20 % load):
   - the source's own connectors (0.03-0.06 N/mm at those stretches): 0.30;
   - the lofted fan that stands in for them (E = 21 MPa, 0.125 mm, roughly 1-2 N/mm in its plane):
     0.57 (60 % load). The chains then carried up to 1.2 N of compression without buckling.

   **The settle has to stay dynamic when the model contains a soft membrane.** A two-step input (dynamic ramp
   + settle to t = 1.2, then STATIC with the loads held) diverged at the first static iteration
   (displacement norm ~2900, 679 negative jacobians): without inertia the wrinkled soft loft's near-zero-
   stiffness modes are singular. Restarts help less than they seem. `febio4 -i m.feb -dump=1 m.dmp` writes
   the last converged state (~95 MB here, overwritten each step), and `febio4 -r m.dmp` (or a
   `<febio_restart version="2.0"><Archive>m.dmp</Archive>...` file) resumes it. The restart file can
   redefine tolerances, `max_refs` and load curves (e.g. the damping curve), but not `time_steps`. An
   added `<Step type="solid"><Control><analysis type="static"/>...` (attribute form; `<analysis>STATIC
   </analysis>` is rejected there) runs only after the original step has ended. So a restart can't cut
   a long settle short. Use a two-step input file instead (`<Step><step id="1">` + `<step id="2">`,
   top-level `<Control>` removed; verified in 4.13). Speed: an iteration of the 54k-equation model
   costs ~1 s. 16 threads were only 1.4x faster than 4. Most of the time goes into failed attempts
   (1000+ with the source's heavy truss mass, 58 with that mass x10), so the path matters more than the
   threads.

   Bending stiffness is the other route: a zig-zag of spring length L under compression F needs roughly
   EI >= F L^2 / 4 (about 3.5 N mm^2 here). A shell tube lofted around the chain cannot add bending on
   its own terms: EI/EA = R^2/2, so on R = 1 mm any useful EI also adds tens of percent of axial
   stiffness in parallel with the springs. If the tube replaces the springs instead, it needs the
   source's axial law. The reference model's old arcus tube had an Ogden fit giving about 1/5 of the
   source truss stress. FEBio 4.13 has elastic beams (`linear-beam`, `febio-xml-format.md` gotcha 25),
   where A and I are set separately. Because the 2-node element shear-locks, a two-element zig-zag is
   resisted by G*A_s, so that must exceed the compression. In the reference model, though, a beam in
   parallel with the tied chain made every static variant fail earlier (t = 0.02-0.13 against 0.29
   without it), and the divergence started away from the arcus. Diagnose that interaction in a mini
   model before counting on a beam.

   Results of the direction check, one change each on the fan-supported run (0.57): the diagnostic's
   failed attempts ran away at a chain node whose drift was almost entirely *out of the fan's plane*.
   A membrane fan supports the chain only in its own plane, so the next weakest direction takes over.
   The same happened dynamically (chain mass x10) at t = 0.75.

4. **Report a comparison table, not just the final answer**, when handing
   this off or writing it up -- "steps completed / time reached / failure
   type" for every variant tried is what lets someone else (or a future
   session) avoid re-deriving conclusions that already exist.

## When to stop guessing and ask for human input

If several structurally different techniques (contact isolation, ramping,
solver settings, continuation methods, stiffness scaling) all fail to
improve on a known baseline, or a "fix" for one failure mode produces an
even more severe failure (e.g. a modest failure becomes a NaN blowup), that
is a signal to stop iterating blindly and either (a) report the specific,
characterized failure and the options tried, or (b) get eyes on the actual
deformed shape at the failure point in a GUI, which nothing in this
text-only workflow can substitute for. Visually inspecting where and how
the mesh is deforming right before failure is often the fastest way to
identify a specific local cause (a self-intersection, a single bad element,
a wrongly-modeled joint) that a log-file-only diagnostic loop cannot see
directly.
