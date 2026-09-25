# Abaqus .inp format reference (as source-of-truth for a FEBio conversion)

`.inp` is a keyword-and-data-line text format, not XML -- much easier to
`grep`/`sed`/`awk` than FEBio's XML, but has its own structural conventions
that matter when using it as the ground truth to validate or correct a
converted FEBio model.

## Structural layout

```
*Heading
*Preprint, ...
** PARTS
*Part, name=PartA
  *Node                          <-- node IDs are PART-LOCAL here
    1, x, y, z
  *Element, type=T3D2            <-- 2-node truss/wire
    1, 1, 2
  *Nset, nset=SetName, internal, generate
    <start>, <end>, <step>       <-- inclusive range, not a literal list
  *Elset, elset=SetName, internal
  ** Section: ...
  *Solid Section, elset=SetName, material=MatName
    <cross-sectional area>,      <-- for a truss/wire section, this one
                                     number IS the area
*End Part
*Part, name=PartB
  ...
*End Part
...
*Assembly, name=Assembly
  *Instance, name=PartA-1, part=PartA
  *End Instance                  <-- no translate/rotate lines = zero offset,
                                     instance coords == part-local coords
  *Instance, name=PartB-1, part=PartB
  *End Instance
  *Nset, nset=SetX, internal, instance=PartA-1, generate
    1, 12, 1
  ** ties, rigid bodies, contact pairs reference INSTANCE-qualified sets
  *Tie, name=Constraint-Name, adjust=no, position tolerance=<mm>
    <slave-surface>, <master-surface>   <-- SLAVE (secondary) FIRST, master (main)
                                            second, per the Keywords Reference
                                            (an earlier version of this note had
                                            the order reversed)
  *Surface, type=NODE, name=SomeName
    <nset-name>, 1.
  *Surface, type=ELEMENT, name=SomeName
    <elset-name>, <SPOS|SNEG>
  *Rigid Body, ref node=<nset>, elset=<elset>, tie nset=<nset>
  ** discrete point-to-point nonlinear spring/joint, NOT a *Tie:
  *Element, type=CONN3D2
    <id>, <InstanceA>.<nodeA>, <InstanceB>.<nodeB>
  *Connector Section, elset=<elset-of-that-one-element>, behavior=<BehaviorName>
    Axial,
    "Datum csys-<name>",
  ** ...repeated per discrete connector element...
  ** visualization-only reference geometry, NOT mechanically active:
  *Display Body, instance=<InstanceName-1>
*End Assembly
** connector force-deflection laws live outside the Assembly, referenced by name:
*Connector Behavior, name=<BehaviorName>
*Connector Elasticity, nonlinear, component=1
  <force>, <relative displacement>   <-- FORCE FIRST (dependent variable first,
  ...                                    same pattern as *Uniaxial Test Data; an
                                         earlier version of this note had it
                                         reversed -- see the connector bullet)
** reference points / origins defined directly in the Assembly (no part):
*Node
  5, x, y, z                     <-- assembly-level node 5; referenced as a bare
                                     number ("24, 6, VW-PeB-1.43"), and by any
                                     *Nset in the Assembly WITHOUT instance=
** BOUNDARY CONDITIONS
*Boundary
  <nset>, PINNED                 <-- fixes all active translational DOFs
  <nset>, ENCASTRE                <-- fixes everything including rotations
  <nset>, XSYMM                   <-- symmetry BC: fixes the DOF normal to
                                      the symmetry plane + 2 rotations, NOT
                                      a full pin -- do not treat XSYMM as PINNED
  <nset>, 1, 1                    <-- explicit DOF range form (dof1_start,
                                      dof1_end[, magnitude])
** MATERIALS
*Material, name=MatName
  *Density
    <value>,
  *Hyperelastic, marlow, poisson=<v>
  *Uniaxial Test Data
    <nominal stress>, <nominal strain>    <-- STRESS FIRST, STRAIN SECOND
    ...
```

## Things that matter for a conversion/comparison task

- **Node IDs restart at 1 for every `*Part`.** They only become globally
  unique once placed into the `*Assembly` via `*Instance`, and even then
  they're conventionally referred to as `InstanceName.NodeID`, not a flat
  global number. A merged FEBio file renumbers everything into one flat
  space -- matching an Abaqus part's original node back to its FEBio
  equivalent has to be done **by exact coordinate**, not by ID (see
  `geometry-cross-referencing.md`).

- **`generate` on `*Nset`/`*Elset` is `(start, end, step)`, inclusive.** A
  quick way to miscount a boundary condition's extent is to read the two
  numbers as a literal 2-item list instead of a range.

- **`*Tie` can pair a `type=NODE` surface (bare node cloud, no facets
  needed) with a `type=ELEMENT` surface (real facets), or two `type=NODE`
  surfaces.** This is why a 1D truss/wire part can be tied to a shell
  surface in Abaqus with zero extra geometry on the truss side -- the tie
  only needs its *nodes*. FEBio's tied-contact equivalent has no direct
  node-only surface type (see `febio-xml-format.md` gotcha 5), so this
  exact pattern needs a workaround when converting/re-deriving it.

- **`*Connector` elements (`CONN3D2` etc.) are a completely separate
  mechanism from `*Tie`** -- a discrete point-to-point joint/spring between
  two specific nodes (often on two different instances), not a
  surface-based constraint. Structurally: one `*Element, type=CONN3D2` line
  per connector (`id, InstanceA.nodeA, InstanceB.nodeB`), a
  `*Connector Section, elset=..., behavior=<BehaviorName>` immediately
  after each one naming its force law, and the force law itself defined
  once, elsewhere in the file, as `*Connector Behavior, name=<BehaviorName>`
  followed by `*Connector Elasticity, nonlinear, component=1` and a
  table of `force, relative displacement` rows (**force first**). A quick
  consistency check: behaviours in one family usually share the same
  displacement grid in column 2 (e.g. 4.7, 9.515, ...) while column 1
  scales with each behaviour's stiffness factor, which only makes sense
  if column 1 is force. A source model that attaches a structure to
  its surroundings via a *sparse handful* of these (rather than one
  continuous tied surface) is describing discrete spring-like attachment
  points, not a continuously-bonded interface -- converting that pattern
  into one continuous tied/bonded surface changes the mechanical intent,
  not just the representation (see the note on synthesized "fan"/"tube"
  geometry a few bullets down).

- **`*Display Body` marks a part instance as visualization-only, with no
  mechanical role in the analysis at all** -- no stiffness, no load path, not
  tied to anything, contributes nothing to the solve. A part with this as
  its *only* `** Constraint:` entry (see the next bullet) is reference
  geometry, full stop; don't assume it's "just missing its connections" and
  try to wire it up as if it should be load-bearing.

- **Search `** Constraint:` comment markers to enumerate every constraint on
  a part, not just whichever type you happen to grep for first.** A part
  can have `*Tie`, `*Rigid Body`, `*Coupling`, `*Display Body`, or an MPC
  attached, and these use structurally different keywords -- grepping only
  for `*Tie` (or only for `CONN3D2` connector elements) and finding nothing
  does **not** mean a part is unconnected; it may be tied or rigid-body-
  constrained by a different mechanism entirely. `grep -n -A1 '^\*\* Constraint:' model.inp`
  lists every named constraint in the `*Assembly` block regardless of its
  underlying keyword, in one pass -- do this before concluding any part is
  mechanically isolated.

- **A `*Rigid Body` whose reference node is fixed fixes the whole part.** In
  `*Rigid Body, ref node=_PickedSetA, elset=_PickedSetB, tie nset=_PickedSetC`
  every node of the elset's elements and of the tie nset moves with the ref node.
  `_PickedSetA` is often `*Nset, nset=_PickedSetA, internal` with no `instance=`,
  i.e. assembly-level node numbers (reference points, `*Node` inside the
  `*Assembly`). A `*Boundary` on another `_PickedSet` holding the same assembly
  node, in all 6 DOFs, then fixes every node of the part, and its section
  (a `*Beam Section` with E 2.1e8, say) is irrelevant. In the reference model the
  CL/USL anchor lines were such rigid bodies. The conversion dropped the beams and
  the rigid bodies, left the line nodes (where 82 connectors anchor) as orphans, and
  pinned the lofts' own edge nodes 0-0.5 mm away instead. Modelling those
  connectors as springs meant fixing the line nodes again.

- **An `*Nset` in the `*Assembly` with no `instance=` refers to
  assembly-level nodes**, meaning nodes defined by `*Node` lines placed directly
  inside `*Assembly` (usually reference points: connector anchors, BC
  origins). It does *not* refer to some part's nodes, even when the numbers
  happen to fall inside a part's node range, and even when an adjacent
  declaration with a similar name carries `instance=`. Before assigning such
  a set to anything, list the assembly-level nodes:
  `awk '/^\*Assembly/{a=1} /^\*End Assembly/{a=0} a&&/^\*Node/{getline; print}' model.inp`.
  Connector endpoints written as bare numbers (`24, 6, VW-PeB-1.43`) are
  these assembly nodes too. **Correction of an earlier version of this note**,
  which suggested such sets belong to a neighbouring instance. Acting on that
  pinned 72 levator-ani nodes (`Set-PM-origins` 7-32 etc. read as LA part
  nodes) when the sets were really 79 assembly reference points. The
  correct LA pins came from the instance-qualified sets
  (`*Nset, ..., instance=LA-new-1`), and those matched the source exactly by
  coordinate.

- **Resolve every `*Tie` into actual node pairs before converting it.** The
  data line is `slave, master`. With a node-based master surface
  (`*Surface, type=NODE`, often a CAE `_CNS_` surface built from a truss's
  nodes) there are no master facets, so each slave node within the
  `position tolerance` is tied to its nearest master node (this is how the
  pairs were resolved in the reference session; if exact Abaqus semantics
  matter, check with a tiny Abaqus model). Resolve the slave side from the
  *actual* surface definition (an element surface's elset can span much more
  of a part than a converter assumed). In the reference case the slave
  surface was all 1170 elements of the muscle, not just one region of it,
  and 82-84 LA nodes per side tied along the arcus instead of the handful a
  one-region FEBio contact could reach. Reproduce node-to-node ties in FEBio with
  linear constraints (`febio-xml-format.md` gotcha 18), not with a
  proximity contact that snaps gaps shut (gotcha 17).

- **Check the `*Step` before assuming the source is a static analysis.**
  `*Dynamic, Explicit` (often a 1 s step with `*Bulk Viscosity` and a
  `SMOOTH STEP` amplitude) is a quasi-static *explicit* analysis. It never forms a
  stiffness matrix, so it tolerates things a static implicit FEBio solve
  cannot: floppy membranes, cable/truss chains with free interior nodes,
  dangling flaps and near-mechanisms. When a converted static FEBio model
  struggles with exactly those parts, that is why. Note that FEBio implicit
  dynamics with the source densities was *worse* in the reference session;
  the fix was removing the near-mechanisms and pinched elements. Densities in
  explicit models are often mass-scaled on some parts (e.g. 1e5x on thin trusses),
  so don't copy them into a FEBio dynamic run unexamined. **Also grep the
  step for `*Contact Pair`.** Contact pairs live inside the step, not the
  assembly, so they are easy to miss in a conversion (reference case: a
  frictionless PVW-to-levator-ani contact existed only there).

- **Count every connector that touches the part you are converting.**
  `awk '/^\*Element, type=CONN3D2/{getline; print}' model.inp` lists every
  connector's endpoints. Tally them per instance and per behaviour. In the
  reference case the FEBio conversion represented only 16 of 40 LA-to-perineal-body connectors (the other 24
  were a "fan" left unattached). The endpoints existed exactly in the FEBio
  mesh (0.0000 mm), so all 40 could be modelled directly as FEBio
  `nonlinear spring`s with points `(elongation, force)`.

- **Connector force-law extrapolation.** Abaqus `*Connector Behavior` has an
  `EXTRAPOLATION` parameter (CONSTANT by default): outside the table the
  force is held constant, so with a table starting at `0, 0` there is no
  compression force at all. FEBio's `<extend>extrapolate</extend>` instead
  continues the first/last segment linearly (compression resistance, and a
  growing force beyond the table). Choose deliberately, and say which you used.

- **`*Boundary` keyword shortcuts are not interchangeable.** `PINNED` and
  `ENCASTRE` are effectively full fixations (verify which DOFs are
  "active" for the element type in question), but `XSYMM`/`YSYMM`/`ZSYMM`
  are **partial** -- a translation-normal-to-plane + two rotations, leaving
  the other two translations and one rotation free. Converting an `XSYMM`
  BC into a FEBio "zero displacement" bc with all three translational DOFs
  fixed silently over-constrains the model relative to the original design
  intent. (Shell formulations that don't expose explicit rotational DOFs
  can only approximate the translational part of the symmetry condition --
  document that as a known simplification rather than silently treating it
  as equivalent.)

- **`*Hyperelastic, marlow` + `*Uniaxial Test Data` on a *continuum* (shell
  or solid) section** is a full 3D strain-energy fit from 1D test data.
  Marlow's deviatoric energy depends on **I1 only**, so refit it to an
  I1-only law: FEBio `Yeoh`, with `scripts/fit_yeoh.py`. A uniaxial-only
  Ogden fit can match the uniaxial curve and still be 2-3x too soft in
  biaxial stretch. With `poisson=`, Marlow also keeps nu constant at every
  strain, which a single FEBio `k` does not. Both points are in
  `febio-xml-format.md` gotcha 23. Otherwise, if the target solver doesn't
  support Marlow directly, this needs
  re-fitting to an available analytical hyperelastic law (Ogden, etc.)
  using the same data points, which is a genuine fidelity approximation
  worth flagging, not a free conversion. `scripts/fit_ogden.pl` in this
  skill does this fit -- an N-term (default tries 1 through 3, prints all so
  you can compare) Ogden least-squares fit via multi-start coordinate
  descent, no external numerical library needed. It reports R^2, a
  Drucker-stability flag (`mu_i*alpha_i > 0` for every term) and a
  monotonicity flag (does the fitted curve ever dip or reverse over the
  data's own strain range) for each N -- watch the monotonicity flag
  especially: it can catch overfitting a plain R^2 comparison misses,
  particularly once N reaches 3 terms (6 free parameters) against a
  data table with only a handful of points. A concave-down/softening
  stress-strain shape (large initial slope, flattening out) commonly has NO
  fully stable fit at any N -- that's expected, not a bug, and worth
  documenting as an extrapolation caveat rather than forcing an artificially
  restricted (positive-only) search to chase stability at the cost of
  underfitting. See `febio-xml-format.md` gotcha 16 for converting the
  resulting `(mu_i, alpha_i)` pairs into FEBio's `(c_i, m_i, k)` convention,
  including the bulk-modulus/Poisson's-ratio relation for `k`.

- **The exact same material definition on a *truss/wire* section**
  mathematically reduces to a plain nominal-stress-vs-nominal-strain
  force law -- no 3D fit needed. Convert it directly: `force = stress x
  cross_sectional_area` (the area is the single number under
  `*Solid Section` for a truss), `x-axis = strain` (dimensionless, so one
  shared table works for every segment regardless of individual length --
  see `febio-xml-format.md` gotcha 7 on FEBio's `measure="strain"` option
  for exactly this case).

- **Watch for parts that look identical in purpose but got separately
  meshed/refined during conversion** -- e.g. a slender ligament-like
  structure modeled as a bare 1D truss chain in Abaqus can get converted
  to a hollow/ribbon "tube" shell surface in a downstream tool so that it
  has real facets to serve as a contact/tie surface. That conversion
  necessarily displaces the surface nodes radially off the true
  centerline by roughly the tube's assumed radius/thickness -- which can
  quietly inflate the gap that a tied-contact tolerance has to bridge.
  Cross-check: compute the actual offset distance between the truss's
  original centerline nodes and the converted surface's nodes; if it's on
  the same order as a contact's `position tolerance`/`max_distance`, that
  offset is probably a contributing cause of contact-engagement
  instability downstream (see `convergence-debugging.md`).

- **A converted geometric construct may not correspond to any single
  Abaqus part at all** -- it can be a downstream tool's *synthesized
  approximation* of a sparse pattern from the source model, most commonly a
  continuous "fan"-shaped patch standing in for what was originally a
  handful of discrete `*Connector` attachment points (previous bullet's
  tube-from-truss case is the same phenomenon for a 1D chain). Two
  consequences worth checking for before assuming a converted domain is a
  faithful 1:1 translation: (1) its node count and boundary shape won't
  match anything in the `.inp` directly -- cross-reference by tracing which
  *original* connectors/nodes its boundary was built to approximate,
  usually by coordinate match against the connector endpoints, not by
  looking for a same-named `*Part`; (2) if the source model attached that
  region via *N* independent connectors, each with its own (possibly
  different) stiffness curve, and the conversion merged them into one
  continuously-tied surface, that is a real, deliberate simplification of
  the mechanical intent (independent point springs -> one bonded interface),
  worth flagging explicitly rather than treating as a value-neutral format
  change -- see the fidelity-approximation framing in `SKILL.md` stage 3.

  **In the reference project these were deliberate modelling choices, not
  converter output. Call them lofts.** The modeller lofted a shell surface over
  the Abaqus line elements (the `CONN3D2` connector lines, trusses or springs
  between a ligament/arcus line and the tissue it holds), then meshed it with
  MMG (the loft tool is not recorded). The model's domain names end in `_fan`
  (e.g. `P-arcus-L_fan`), but the modeller calls these surfaces lofts and is not
  sure where "fan" came from, so use "loft" in writing. The "tubes" were made
  separately, as connection surfaces: the anchor lines that fixed BCs act on,
  or the lines the LA ties to (e.g. R = 1 mm, t = 1/(2 pi), so the wall area
  equals a 1 mm^2 truss section). Most tubes were later removed (their nodes
  fixed, or replaced by truss chains), and the modeller is wary of them
  because they change the response a lot. Treat a loft as the intended
  representation of its connector family, and do three things:
  1. **Check it is still anchored at both ends.** When the structure a loft was
     lofted onto is replaced (tube -> truss chain), the loft's edge can be left
     free. The loft still shares nodes with the tissue, so it looks connected.
     In the reference model the posterior-arcus lofts hung loose from the
     vaginal wall through three model versions, so the 26 source connectors
     `P-arcus-L/R-1..13` had no counterpart at all. For each version, list
     what holds every boundary node of each loft (other domains, BCs, springs,
     contact surfaces, linear constraints) and compare across versions
     (`claude_diag`-style `parcus_fan_history.py`). Reattach it by sharing
     nodes with the line it was lofted over, not by hidden linear constraints:
     FEBio Studio doesn't draw constraints, so a constrained loft still looks
     loose. Then check *every* loft edge node along that line, not just the ones
     at the line's own nodes. After the merge, 8 loft edge nodes per side still
     lay on the arcus chain line without being attached. Under load they
     scalloped 5-10 mm away from the chain, and the gaps showed in FEBio Studio.
     If the line is a chain of strain-measure springs, split the spring through
     each such node, which keeps the law. The same goes for a loft edge pinned
     only at the source's point anchors (a BC on the connector origins): pin every
     loft node along that edge (they lie on the straight segments between the
     anchors). On a weak fitted loft the free ones swung 3.6 mm. Pinning them gave
     the same solution with half the failed attempts, in two runs. Also check the
     tied lofts' ties actually hold (`febio-xml-format.md` gotcha 17: `penalty` is
     per node).
  2. **Give every loft its own material, fitted to the source elements it
     replaces.** Don't reuse a generic or shared material. In the reference
     model all lofts shared one isotropic-elastic E = 21 MPa, 0.125-0.49 mm
     material. For the posterior-arcus loft that was 16-60x stiffer than its 13
     soft connectors (0.117 N at 4.7 mm, 5.9 N at 47 mm, a J-shaped curve).
     Reattached, it halved the LA displacement at the same load, dragged the
     vaginal wall 17 mm, and stopped the dynamic run at 84 % load, where the
     connector version finished. Fitting method (strip model): connector i
     (length L_i, from the source nodes) becomes a loft strip of width w_i (its
     spacing along the two loft edges) and the loft thickness t. Pulling the
     anchored edge a distance u away gives stretch 1 + u/L_i. Fit a hyperelastic
     law so that sum_i P(1 + u/L_i) w_i t matches sum_i F_i(u) over the working
     range of u. A 1-term Ogden (FEBio c1 = 2 mu, m1 = alpha) follows a J-shaped
     connector curve well: c1 = 0.259 MPa, m1 = 3.5 matched the 13 connectors to
     2.8 % rms over u = 2-40 mm (`claude_diag`-style `fit_fan_to_connectors.py`).
     Linear or neo-Hookean laws cannot stiffen enough. A continuous sheet is
     somewhat stiffer than independent strips, so check the loft's actual pull in
     a run. **Keep the strips; don't collapse them into one cross-section.** A
     fit of the total connector force over the loft's average cross-section (area /
     mean length x t), at strain u / mean length, came out 1.5-2x stiffer where the
     model worked, and in situ it moved the ends further from the connectors (USL
     60 % of the connectors' opening vs 67-79 % for the strip fit). The real loft
     strains its short part more (u / its own length), so with a stiffening law the
     short end pulls harder, and only the strip model builds that in. The
     cross-section curve is still a readable "tissue stress-strain" to report. Two
     more settings decided whether the loft ran:
     - **k from the stiffened modulus, not mu0.** With k = 50 mu0 the loft swelled
       (J median 1.32, max 4.25) once the m1 = 3.5 law had stiffened ~5x at ~3x
       stretch, and the run failed at t = 0.96. With k = 250 mu0 J stayed <= 1.5 and
       it finished.
     - **The thickness of the model's other lofts** (0.49 mm), with c1, k and density
       scaled by t_old/t_new to keep the fitted pull and the mass. A 0.125 mm soft
       sheet in a tension field (stretched ~3x, squeezed crosswise) wrinkles and
       crawls.
     Result: the fitted loft reproduced the connector run like for like at full load
     (LA median / p90 / max 29.5 / 40.3 / 42.8 mm vs 29.7 / 40.1 / 42.8).
  3. **Compare with the exact alternative:** the connectors themselves as FEBio
     nonlinear springs (force-first table swapped, `extend` constant =
     tension-only). Run both, one change each. The sharpest in-situ check is the
     connector-end elongation (distance change between each connector's two end
     nodes) in the loft run against the spring run at the same times. The
     posterior-arcus loft at c1 = 0.066 gave median 17.3-18.0 / max 25.8-26.6 mm
     at t = 1, against 18.0-18.5 / 24.8 mm for the springs. The LA displacement was
     no test, because the link moved it by only 5-10 %.
  4. **Audit every loft against its own connectors, not just the one in trouble**
     (`claude_diag`-style `loft_survey.py`, then `loft_fit_all.py`). Map every
     `CONN3D2` family by both ends (0 mm coordinate mismatch), say which loft
     carries it, and list how each loft is held: shared nodes, `tied-node-on-facet`,
     BC-fixed edge nodes, or a spring that already carries a connector exactly
     (don't fit that one again). Widths for the strip model: order the connectors
     along the tissue edge. Consecutive connectors bound a ruled quad; each owns half
     of each neighbouring quad; scale so that sum w_i L_i = the loft's mesh area.
     Spacing-based widths counted the posterior-arcus loft area twice (676 vs
     325 mm^2). That tool's c1 still matched in situ, where the area-consistent fit
     was 1.4x stiffer: a continuous sheet is stiffer than independent strips
     (planar tension explains only 1.0-1.2x of that). What the reference audit found:
     - the other 12 lofts all used one isotropic-elastic E = 21 material: the Abaqus
       CL/USL *pipe beam* material (`*Beam Section, section=PIPE`, r 0.5, wall 0.49
       = the lofts' thickness; Abaqus E 2.1e8, a rigid anchor line) scaled down.
       With the beam's density 7.8e-07 the lofts weighed 1.65 kg against 69 g for
       all the tissue, which matters in dynamic runs;
     - against their connectors: AVW-Para 2.5x, CL 40-60x, USL 180-350x too stiff;
       PM 0.6-0.8x (too soft); five lofts stood in for connectors impaired to almost
       nothing (6e-05 N at 4.7 mm, 0.003 N at 47 mm), so they were 10^4-10^7x too
       stiff, and one was even carrying compression, which the source connectors can't;
     - one material per loft can't share load the way the connectors do. Every
       connector in a group has the same force-elongation whatever its length, but a
       strip's force follows u/L. Per-group fits differed 4-14x within a loft;
     - a 1-term Ogden can't follow a long slack toe (0.2 N at 25 mm; rms 22-24 %);
     - a loft fitted to near-zero connectors (E0 ~ 60 Pa here) must also get a light
       density (tissue, or near zero; connectors are massless). At the inherited
       density one such loft weighed 64 g (all the tissue: 69 g), and in a dynamic
       run its inertia would dwarf the connector forces (0.0005-0.004 N). It also
       stretches hugely (strips of short connectors ~4x at 5 mm, ~12x at 20 mm), so
       check its J for swelling (in the runs it stayed 0.95-1.15: no swelling).
  5. **Then check all the fitted lofts in situ at once:** one run with the fitted
     lofts, one with all their connectors as springs (`claude_diag`-style
     `variants6.lofts_to_connectors`: remove the loft and any tie on it, add the
     springs, re-fix orphaned anchors), same base, connector-end elongations at the
     same time. What the reference check found (175 connectors, 12 families):
     - the near-zero families matched their connectors (PM_PeB median / max
       2.4-2.6 / 6.6 vs 2.5-2.7 / 6.6 mm), and so did the posterior arcus;
     - CL medians matched with the max at 70 %; USL and a small 2-connector loft
       reached 55-70 %, so a loft can be ~1.5-2x stiffer in situ than its strip fit;
     - **a loft carries compression and the connectors don't** (a table from (0, 0)
       with constant extrapolation is tension-only). Where the tissue moves toward the
       anchor, the connectors shorten freely (AVW-Para: median -4.7 mm) while the
       fitted loft resists (-1.1; the E 21 loft -0.2). No material fit can fix that.
       Look at the sign of the elongations in the spring run before fitting;
     - the E 21 lofts had held every end within ~1 mm. Faithful supports let the
       apex open 3-9 mm, yet the LA moved only +-4 %. Find out what the output of
       interest is sensitive to before tuning lofts further;
     - the spring version converged best (0 negative Jacobians with the line
       search), so where the source has connectors, springs are a usable model in
       their own right. Lofts fitted to near-zero connectors made convergence worse
       (a membrane of almost no stiffness, with hundreds of interior nodes) and
       passed only with the line search on.

- **Audit connector families by both ends, not by the part under study.** A
  connector audit that asks "is any connector to part X missing?" misses a
  family attached to what X hangs from. List every `CONN3D2` family with its two
  end parts, then name the FEBio construct that carries each family in the
  *current* version (springs, a fan, a tie) and check that it is attached at
  both ends. In the reference session an LA-scoped audit ("only the 40
  sphincter connectors touch the LA") was correct but missed the 26 arcus-to-
  vaginal-wall connectors on the line the LA hangs from.

- **Cross-check load magnitudes, not just BCs/contacts/materials.** `*Dsload`
  (distributed surface load -- pressure/traction) is an easy thing to leave
  out of a conversion-fidelity pass, since it's a single number per surface
  and doesn't have the structural complexity that draws attention the way a
  BC or contact does. Confirmed real case: an Abaqus source had two surfaces
  loaded with the *same* pressure magnitude (`*Dsload, amplitude=Amp-1` on
  both, one value shared), and the converted FEBio file's `<surface_load>`
  block for one of them had silently gone to `0` -- passing every other
  fidelity check while quietly changing the loading scenario. Diff the
  `*Dsload` lines (surface name, magnitude, amplitude reference) against
  FEBio's `<Loads>` block's `<pressure lc="...">` values directly, the same
  way BCs and materials get cross-checked, rather than assuming a load that
  "looks present" on both sides also matches in magnitude.

- **A `*Boundary` BC can survive a conversion as an unreferenced `<NodeSet>`.**
  The node group exists in the FEBio file under a matching name, with exactly
  the right members (verified by coordinate), but no `<bc>` in `<Boundary>`
  ever uses it -- a different, more silent failure mode than an orphaned BC
  whose `node_set` points at nodes a remesh moved (`geometry-cross-referencing.md`
  Rule 12, `scripts/bc_audit.py`): here the constraint itself was simply never
  written. The symptom is not an obvious defect. Confirmed case: an Abaqus
  `XSYMM` on 328 vaginal-wall midline nodes converted to a correctly-named,
  correctly-populated FEBio `<NodeSet>` that nothing referenced; the nodes
  were instead held only by unrelated low-stiffness ground springs. The model
  ran, converged cleanly to 60 % of the intended load with small,
  physically-plausible-looking displacements, then failed -- indistinguishable
  from a genuine physical/mesh limit until every `*Boundary` set name was
  cross-checked against FEBio's own `<Boundary>` block, not just against its
  `<NodeSet>` declarations. Restoring the one missing `<bc>` took that model
  from 60 % of load to 97 %+.

  Then check what the stand-in left behind. Those ground springs were
  **zero-length** (anchor on the node), so each one resists all three directions,
  unlike the XSYMM (u_x only on solids). With the BC restored and the springs
  kept, the midline moved 14 % less in-plane than its neighbours: a visible
  groove. The springs also carried ~20 % of the net applied load. "The symmetry
  plane looks stiff in y and z" was the springs, not the BC. Nodes beside the
  plane moved only 0.2 mm sideways, so the BC itself was doing little.

  ```bash
  # Abaqus *Boundary nset names
  awk '/^\*Boundary/{getline; print $1}' model.inp | tr -d ',' | sort -u
  # FEBio node_set attributes actually used by a <bc>
  grep -o 'node_set="[^"]*"' model.feb | sort -u
  # any Abaqus name missing from the second list has a NodeSet but no <bc>
  ```

## Extraction commands that hold up on huge `.inp` files

```bash
# List every part name
grep -n '^\*Part, name=' model.inp

# Pull one part's full node+element+section block (handles quoted names)
sed -n '/\*Part, name=PartName/,/\*End Part/p' model.inp

# Pull one material's full definition
sed -n '/\*Material, name=MatName/,/\*Material,/p' model.inp | head -n -1

# Find every Tie constraint and what it connects
grep -n -A2 '^\*Tie,' model.inp

# Resolve which instance a *Surface/*Nset actually belongs to
grep -n 'nset=_PickedSetNNN,' model.inp   # look for ", instance=" on that line

# List EVERY constraint on any part, regardless of keyword (Tie, Rigid Body,
# Coupling, Display Body, MPC, ...) -- the only reliable way to confirm a
# part is genuinely unconnected rather than just missing from a *Tie grep
grep -n -A1 '^\*\* Constraint:' model.inp

# Find every discrete *Connector element and what it joins (distinct from *Tie):
# endpoints as "Instance.node" or a bare number (= assembly-level reference node)
awk '/^\*Element, type=CONN3D2/{getline; print NR": "$0}' model.inp
# ...and tally endpoints per instance
awk '/^\*Element, type=CONN3D2/{getline; print}' model.inp | awk -F, '{print $2; print $3}' | sed 's/\..*//' | sort | uniq -c
grep -n '^\*Connector Behavior, name=' model.inp

# Assembly-level reference nodes (what instance-less *Nsets and bare-number
# connector endpoints refer to)
awk '/^\*Assembly/{a=1} /^\*End Assembly/{a=0} a&&/^\*Node/{getline; print}' model.inp

# Analysis procedure and contact pairs (both live in the *Step, not the Assembly)
awk '/^\*Step/{f=1} f&&/^\*(Static|Dynamic|Contact Pair|Dsload|Amplitude)/' model.inp
```
