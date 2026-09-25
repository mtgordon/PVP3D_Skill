# FEBio XML (.feb) format reference

FEBio's `.feb` format (febio_spec version 4.0 in the reference session; check the
`<febio_spec version="...">` root attribute, since parameter names shift between
major versions) is XML, but its parser has several undocumented strictness rules
that only surface as a "Reading file ... FAILED!" error with a line number. This
file collects everything learned the hard way. Read the top-level section layout
first, then the gotcha list before editing any real model.

## Top-level section layout

```
<febio_spec version="4.0">
  <Module type="solid"/>
  <Control>        analysis type, time stepper, nonlinear solver settings
  <Globals>        gravity, referential constants
  <Material>       <material id=".." name=".." type="..">...</material> blocks
  <Mesh>           ALL geometry: Nodes, Elements, NodeSet, Surface, SurfacePair,
                   DiscreteSet   <-- yes, DiscreteSet lives here, not in <Discrete>
  <MeshDomains>    <SolidDomain>/<ShellDomain name=".." mat="..">, binds an
                   Elements block to a material
  <Boundary>       <bc name=".." node_set=".." type="..">...</bc>
  <Loads>          <surface_load>, pressure/traction, referencing a load curve
                   via lc="<id>"
  <Contact>        <contact name=".." surface_pair=".." type="..">...</contact>
  <Constraints>    <constraint type="linear constraint"> exact node-to-node ties
                   (gotcha 18) -- NOT inside <Boundary>
  <Discrete>       <discrete_material id=".." name=".." type="..">...</discrete_material>
                   blocks, THEN <discrete dmat="<id>" discrete_set="<name>"/>
                   binding lines (mirrors how MeshDomains binds Elements to
                   Material -- discrete_material is to <discrete> as material
                   is to ShellDomain)
  <LoadData>       <load_controller id=".." name=".." type="loadcurve"> curves
  <Output>         plot/log file variable requests
</febio_spec>
```

Everything is referenced **by name**, not by position: `node_set="X"` looks up a
`<NodeSet name="X">`, `surface_pair="X"` looks up a `<SurfacePair name="X">`
which itself has `<primary>`/`<secondary>` pointing at `<Surface name="...">`
blocks, `discrete_set="X"` looks up a `<DiscreteSet name="X">`, `mat="X"` looks
up a `<material name="X">`. This makes cross-file renames easy to get wrong
silently -- a typo just means the referenced block doesn't get used, with no
error.

## Gotchas (in order of how much time each one cost)

### 1. Node IDs must increase monotonically across the ENTIRE `<Mesh>` section

This is the single most expensive gotcha in this domain. It is **not** just a
per-`<Nodes>`-block rule -- it is global across every `<Nodes>` block in the
file, regardless of what other tag types (`<Elements>`, `<DiscreteSet>`, etc.)
appear between them.

Symptom: `tag "node" (line N) : invalid value for attribute "id"` -- even though
the actual `<node id="...">...</node>` line at N is syntactically perfect and
is not a duplicate ID.

Diagnosis: check whether the id on the failing line is **lower** than the
highest id declared in any earlier `<Nodes>` block in the file. If new nodes
are inserted in the *middle* of the mesh (e.g. replacing an old part's
geometry in place) but use ID numbers higher than what currently exists,
everything is fine locally -- until the file continues afterward into
pre-existing `<Nodes>` blocks with lower IDs, and *that* triggers the error,
often hundreds or thousands of lines after the actual new content.

Fix: always append brand-new nodes as new `<Nodes>` blocks placed **after**
whatever `<Nodes>` block currently has the highest node ID in the file -- i.e.
at the very end of the node declarations, right before whatever comes next
(usually an `<Elements>` block). Do not insert new high-ID nodes in the middle
of the file just because that's where the *logically related* old content
used to be. Geometry that is logically grouped can still be declared in a
different location in the file than where you conceptually think of it.

This also means: if you build any new content in stages, you may need to move
it *twice* -- once you discover the constraint, both the new `<Nodes>` blocks
**and** anything that references them (`<DiscreteSet>` delem entries,
`<Surface>` facets) must physically sit after all their referenced nodes are
declared, i.e. after your new `<Nodes>` blocks, not before.

Practical check before running FEBio: extract every `<node id="N">` in file
order and confirm N is non-decreasing (see `scripts/check_monotonic_node_ids.sh`
in this skill for a ready-to-use version):

```bash
grep -oE '<node id="[0-9]+">' model.feb | grep -oE '[0-9]+' | awk 'NR>1 && $1<prev{print "DROP at line "NR": "prev" -> "$1} {prev=$1}'
```

### 2. `<DiscreteSet>` lives in `<Mesh>`; the material and the binding live in `<Discrete>`

```xml
<Mesh>
  <DiscreteSet name="my_springs">
    <delem>101,102</delem>
    <delem>102,103</delem>   <!-- multiple delem per DiscreteSet is fine -->
  </DiscreteSet>
</Mesh>
...
<Discrete>
  <discrete_material id="1" name="my_spring_mat" type="nonlinear spring">
    ...
  </discrete_material>
  <discrete dmat="1" discrete_set="my_springs"/>
</Discrete>
```

Putting `<DiscreteSet>` inside `<Discrete>` gives `tag "DiscreteSet" (line N) :
unrecognized tag`. All `discrete_material` definitions come first in
`<Discrete>`, then all `<discrete dmat=".." discrete_set=".."/>` binding lines
after -- this matches the order already used for hundreds of pre-existing
entries in a typical converted model, so grep an existing `<discrete ` line
for the exact pattern before adding your own.

### 3. `<discrete dmat=".." discrete_set="..">` resolution needs a real material
domain to exist somewhere in the file

A minimal test file containing *only* `<Nodes>` + `<DiscreteSet>` +
`<discrete_material>` + `<discrete>` (no `<Material>`/`<MeshDomains>`/
`<Elements>` at all) fails with `tag "discrete" (line N) : invalid value for
attribute "discrete_set"` even though the DiscreteSet is spelled identically
and the name matches exactly. Adding one trivial real material + one
`<Elements>`/`<ShellDomain>` triple anywhere in the file (even unrelated to
the springs) makes the exact same discrete markup parse and run. This looks
like FEBio's mesh/domain initialization needs at least one real domain to run
before it will resolve later cross-references -- irrelevant for any real
converted model (which always has plenty of domains already) but a trap when
building an isolated test file to verify syntax (see gotcha 6).

### 4. `node_set` on a `<bc>` must reference a *named* `<NodeSet>`, never a raw
node ID or a bare list

`<bc name="fix" node_set="1" type="zero displacement">` fails with `invalid
value for attribute "node_set"`. You must first declare
`<NodeSet name="fixset">1</NodeSet>` inside `<Mesh>`, then reference it by
name: `node_set="fixset"`.

### 5. `<Surface>` cannot be built from bare node references -- only facet elements

Abaqus lets a `*Surface` be `type=NODE` (a bare point cloud with no facets) or
`type=ELEMENT`. FEBio's `<Surface>` only accepts facet elements (`<quad4
id="..">n1,n2,n3,n4</quad4>`, `<tri3 id="..">n1,n2,n3</tri3>`, etc.) -- a plain
`<node id="1">5</node>` line inside `<Surface>` gives `tag "node" (line N) :
unrecognized tag`. If the Abaqus source ties a bare 1D structure (a truss/wire
chain with no natural facets) onto a shell surface using a node-type surface,
FEBio has no direct equivalent for the "no facets" side. See
`geometry-cross-referencing.md` for the pattern used to work around this
(a facet built from the truss's own nodes plus a tiny offset "phantom" node
per position, giving a valid near-zero-width ribbon of `quad4` facets to
stand in as the secondary contact surface) -- **and** the follow-on trap that
pattern creates (gotcha 8). **Prefer gotcha 18 instead:** when the Abaqus tie
is node-to-node (or node-to-truss), FEBio linear constraints reproduce it
exactly with no facets, no phantom nodes and no contact at all.

However: a `<Surface>`'s facet elements do **not** need to be part of any real
material `<Elements>`/domain -- a facet built purely for contact/tie geometry,
referencing nodes that carry no material stiffness of their own, parses and
solves fine as long as gotcha 3's "at least one real domain somewhere" rule is
satisfied.

### 6. Verify unknown syntax in a throwaway minimal file, never on the real model

Every gotcha above was discovered by writing a 20-40 line standalone `.feb`
file, running `febio4.exe -i test.feb`, and reading the exact error -- never
by guessing against a multi-MB real file where a failure could be anywhere
and re-running takes much longer to even reach a parse stage. Keep a minimal
skeleton on hand:

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Control>
    <analysis>STATIC</analysis>
    <time_steps>2</time_steps>
    <step_size>0.5</step_size>
    <solver type="solid"><max_refs>15</max_refs></solver>
  </Control>
  <Material>
    <material id="1" name="mat1" type="neo-Hookean">
      <density>1</density><E>1</E><v>0.3</v>
    </material>
  </Material>
  <Mesh>
    <!-- put the construct you're testing here -->
  </Mesh>
  <MeshDomains>
    <!-- a real ShellDomain/SolidDomain if the test needs one (gotcha 3) -->
  </MeshDomains>
</febio_spec>
```

Iterate on this in seconds instead of minutes, and only copy the verified
pattern into the real file once it parses (and, ideally, once it runs to
`NORMAL TERMINATION` for at least a couple of trivial time steps).

### 7. Discovering which element/material types actually exist, with no internet
and no PDF-text-extraction available

If the bundled FEBio manual PDF can't be read (`pdftoppm` missing is a common
failure) and there's no network access, the installed `.dll`s are still full
of the exact XML type-name strings the parser matches against -- because
those strings are compiled-in literals used for factory registration. Extract
them with PowerShell (works even though the files are binary):

```powershell
$f = "C:\Program Files\FEBioStudio\bin\febiomech.dll"   # try febiomech.dll,
                                                          # fecore.dll, febioxml.dll
$bytes = [System.IO.File]::ReadAllBytes($f)
$text  = [System.Text.Encoding]::ASCII.GetString($bytes)
[regex]::Matches($text, "[A-Za-z][A-Za-z0-9_ -]{2,40}") |
  ForEach-Object { $_.Value.Trim() } | Sort-Object -Unique |
  Where-Object { $_ -match "(?i)truss|spring|hyperelastic|whatever_keyword" }
```

This is how it was confirmed that FEBio ships only a `"linear truss"` material
(no nonlinear/tabulated truss material -- see `geometry-cross-referencing.md`
for why that mattered) and that the nonlinear spring's `<measure>` element
accepts `"strain"` in addition to `"elongation"`/`"stretch"` (letting one
shared force-vs-strain table serve every spring in a chain regardless of each
segment's individual length, instead of needing a separately-rescaled table
per segment).

**When you need exact parameter *values*/structure, not just type names, and
dll-string extraction alone won't give you that**: a compiled-in factory
string confirms a type exists and often turns up its bare parameter names as
neighboring strings, but not their expected value ranges, defaults, or exactly
how a `<Surface>`/`<SurfacePair>` needs to be wired up for it. If the user (or
anyone) has a working FEBioStudio project using the construct in question,
ask them to do **File -> Export -> FEBio** and hand over the resulting `.feb` --
a real, solver-validated example is unambiguous ground truth for exact syntax
in a way that no amount of binary-string archaeology or plausible-looking
guessing can match. This is a legitimate escalation step, not a last resort to
feel bad about reaching for: it was faster than continuing to guess, and
resolved a case where the dll strings alone had already led to a wrong
assumption about which side of a contact needed real facets (gotcha 11).
FEBioStudio's own project format (`.fs2` and similar) is a proprietary binary
container, not readable directly -- don't try to parse it; get the plain-XML
`.feb` export instead.

A DLL string can confirm a type exists without meaning it's a good fit: `dll`
strings confirmed FEBio ships a `udg-hex` (uniform-deformation-gradient, 1-point)
solid formulation, which looked like a natural stand-in for an Abaqus `C3D8R`
reduced-integration hex. Tested head-to-head on the same vaginal-tissue Ogden
material in a small compression test, it was **12.6x stiffer** than FEBio's
default (three-field) hex -- not equivalent at all -- and its XML has no
`<hourglass>` tag exposed despite `UDGHourglassForces`/`UDGHourglassStiffness`
strings existing in the DLL (hourglass control, if any, isn't user-facing here).
Don't assume a same-sounding element formulation is a faithful translation of
an Abaqus element without checking its actual stiffness response.

### 8. A geometrically "free" node with zero material stiffness, pulled by a
contact/tie, is a nearly-singular DOF

If you build a facet purely for contact/tie geometry (gotcha 5's workaround)
using nodes that carry no material, no spring, and no boundary condition of
their own, and that facet is then used as the *secondary* side of a tied
contact, the only thing resisting that node's motion in the whole model is
the contact's own penalty stiffness. Combined with any strongly nonlinear
element nearby (a rapidly-stiffening hyperelastic spring, for instance), this
tends to produce catastrophic solver blowups (`NAN detected`, residuals in
the 10^10-10^19 range) rather than a normal convergence failure -- and the
failure shows up instantly, at a vanishingly small load fraction, even with a
heavily ramped-in contact penalty. This is diagnosable by isolation testing
(see `convergence-debugging.md`): disable the tie and confirm the rest of the
new geometry is stable on its own; if it is, and re-enabling *any single* one
of several new ties reproduces the same blowup at the same tiny scale, an
unsupported "phantom" DOF used only for tie geometry is the prime suspect.

The fix that was *attempted but not confirmed working* in the reference
session was linking each phantom node back to its real geometric partner with
a small linear spring, so the tie has real stiffness to react against. It did
not resolve the blowup on the first attempt (made it worse), which means
either the stiffness value chosen was wrong, or the true cause was something
else layered on top (possibly a duplicate-tie constraint at the pinned
endpoint, or over-constraint from multiple ribbon facets sharing a node -- see
`lessons-learned.md`). If picking this up again: try isolating further by
testing a single spring segment with a single tied phantom node in the
minimal test-file harness (gotcha 6) before touching the full model, and
check whether the surface's facets double-count shared nodes across adjacent
quads in a way that creates redundant constraint equations.

**Resolved in a later session -- don't rebuild the phantom ribbon.** The
construction is singular by design. A phantom node held by one spring is
supported in one direction only; the other two directions depend entirely on
the tie. If the tie's penalty is ramped in from zero with a load curve (the
usual mitigation), that support is effectively nothing at t~0. The result is
a residual of about 1e9 on the very first step. Any part of the chain beyond
`max_distance` is left with no support at all. And `tied-node-on-facet`
additionally drags engaged nodes onto the surface (gotcha 17). The
replacement that worked, verified to about 1e-6 mm over 199 ties on a full model, is
gotcha 18: tie each Abaqus slave node to its master node with linear
constraints, and delete the phantom nodes, ribbon surfaces and tied contacts.

**A structurally different case where the spring fix *did* work, confirmed**:
the failure mode above is specifically about a node with *no* stiffness
contribution of its own -- a true "phantom" that exists only for tie geometry.
A separate scenario that looks superficially similar but isn't the same
problem: duplicating a boundary node of an *already-continuously-meshed,
materially-real* domain (e.g. splitting a shell surface's shared-node
connection to a solid it touches, so the two can be tied instead of merged)
and linking the duplicate back to its original position with a discrete
linear spring. Here the duplicated node is *not* a phantom -- it's still part
of a real shell domain with its own genuine material stiffness in every
in-plane direction; the spring only needs to react against out-of-plane/
normal pull toward the solid, not carry the node's entire structural role the
way it would for a true zero-stiffness phantom. This version worked cleanly
(stable convergence, no blowup, spring stiffness on the order of the model's
own softest material rather than needing careful tuning) on a real multi-
domain model. The distinguishing question before assuming either outcome
generalizes to a new case: does the node carry any of its own material
stiffness independent of the tie, or is the tie/spring genuinely the *only*
thing resisting its motion? Only the latter is gotcha 8's near-singular-DOF
failure mode.

### 9. "isolated vertices removed" is (usually) harmless

FEBio prints `WARNING: N isolated vertices removed` for every node declared
that no `<Elements>` block references. If a conversion step leaves orphaned
nodes behind (e.g. deleting an old shell domain but not the nodes it used, or
CAD-export leftovers never cleaned up), this warning fires but has zero effect
on the analysis -- the nodes are silently dropped. Don't chase this warning by
itself; only worry if the *count* jumps unexpectedly after an edit that should
not have orphaned anything, since that can also flag a mistake (e.g.
forgetting to update a reference after renaming a node range).

### 10. `WARNING: No contact pairs found for tied interface "X"` means the tie
is silently inert, not an error

If a `tied-node-on-facet` (or similar proximity-based) contact's
`max_distance` is smaller than the real gap between the two surfaces in the
reference configuration, FEBio finds zero pairs, prints this warning, and
proceeds as if the contact didn't exist -- no error, no failure, just missing
physics. This can hide in a working-looking model for a long time. Always
grep the log for this string after any change to contact tolerances or
geometry; see `convergence-debugging.md` for how to measure the real gap and
choose an appropriate tolerance.

### 11. The three tied-contact types are not interchangeable

FEBio ships (at least) three distinct tied-contact formulations, confirmed
by extracting factory-registration strings from `febiomech.dll` (see
gotcha 7's technique) and mapped to their internal classes:
`tied-elastic` -> `FETiedElasticInterface`, `tied-facet-on-facet` ->
`FEFacet2FacetTied`, `tied-node-on-facet` -> `FETiedInterface`. They differ in
ways that matter when picking one:

- **`tied-elastic`** is sensitive to facet winding and exposes
  `flip_primary`/`flip_secondary` parameters to correct it; it also NaNs on
  a true zero-gap (exactly coincident primary/secondary points) -- needs a
  small deliberate offset between the two sides to avoid a singular contact
  normal.
- **`tied-facet-on-facet`** requires an actual, reasonably well-formed
  facet-based `Surface` on *both* sides -- no bare-node secondary is
  possible even in principle (see gotcha 5).
- **`tied-node-on-facet`** is, empirically, insensitive to facet winding
  (no `flip_primary`/`flip_secondary` parameters exist for it at all) and
  its `search_tolerance` is designed to handle near/exact-coincident points
  natively -- a duplicated tie node can sit at zero offset from its target
  without the `tied-elastic` NaN (**but a *non-zero* initial gap gets closed,
  not preserved** -- gotcha 17 -- so it is only safe when engaged nodes
  already sit on the primary surface or `node_reloc=1` is set). This makes it the more forgiving choice
  when winding correctness can't be fully guaranteed on one side (e.g. a
  converted/synthesized surface, see `abaqus-inp-format.md`) or when the
  tie needs to close a genuinely zero-width gap. **It also has no
  `auto_penalty` parameter at all** (confirmed absent from its full
  parameter list: `laugon`, `tolerance`, `penalty`, `minaug`, `maxaug`,
  `search_tolerance`, `offset_shells`, `max_distance`, `special`,
  `node_reloc` -- compare to `sliding-elastic`/`tied-elastic`, which both
  have one). Its `penalty` is therefore *always* a raw, absolute stiffness
  value, never auto-scaled from the surrounding materials -- a value copied
  from a different model/example (e.g. `penalty=1` from a stiff generic
  test case) can be wildly wrong for a soft hyperelastic tissue model and
  will crush it on the very first iteration; set it relative to your own
  model's own softest engaged material, not by analogy to an example file
  with different materials. Minimal working example (values are a
  starting point for a typically-scaled soft-tissue model, not universal):

  ```xml
  <SurfacePair name="MyTie">
    <primary>TargetSurface</primary>
    <secondary>SourceSurface</secondary>
  </SurfacePair>
  <contact name="MyTie" surface_pair="MyTie" type="tied-node-on-facet">
    <laugon>PENALTY</laugon>
    <tolerance>0.01</tolerance>
    <penalty>0.01</penalty>
    <minaug>0</minaug>
    <maxaug>10</maxaug>
    <search_tolerance>0.0001</search_tolerance>
    <offset_shells>0</offset_shells>
    <max_distance>0.5</max_distance>
    <special>1</special>
    <node_reloc>0</node_reloc>
  </contact>
  ```

Confirm which type a given `.feb` uses before assuming winding matters (or
doesn't) -- verify empirically per gotcha 6's minimal-file approach rather
than assuming based on this list alone, since the winding-insensitivity of
`tied-node-on-facet` was only established by A/B testing single-threaded
deterministic runs before and after a winding fix and observing no change
in the failure point.

`tied-elastic` does not accept every parameter `sliding-elastic` does, even
though they're close relatives in FEBio's contact-interface family --
copying a working `sliding-elastic` block and just changing `type=` is not
safe. Confirmed invalid on `tied-elastic`: `seg_up` (`tag "seg_up" :
unrecognized tag`). When converting a sliding-type contact block to
`tied-elastic`, trim to a conservative core (`laugon`, `tolerance`,
`gaptol`, `penalty`, `auto_penalty`, `two_pass`, `search_tol`,
`symmetric_stiffness`, `minaug`, `maxaug`) and add back only what's
confirmed valid via gotcha 6's minimal-file check, rather than carrying
over a full sliding-elastic parameter block wholesale.

### 12. A proximity-based contact's `max_distance="0"` means *unlimited*, not zero tolerance

This is a specific, expensive-to-discover semantic trap for
`tied-node-on-facet` (and likely other proximity-gated contact types):
`<max_distance>0</max_distance>` does **not** mean "require exact/zero-gap
contact" -- it means *no distance limit at all*. With it left at the default
`0`, every node on the secondary surface ties to whichever primary facet is
nearest, no matter how far away, including nodes that should never have
engaged at all. On a secondary surface with many genuinely-interior/far
nodes and a large primary surface, this can apply large spurious forces
across the whole secondary domain from the very first increment. Always set
an explicit, real distance appropriate to the actual measured gap (see
`convergence-debugging.md`'s point-to-facet distance technique) -- never
leave this parameter at its literal default when non-zero engagement
selectivity is intended.

### 13. A domain touching two different tied-contact targets at a seam can get tied to both at once ("double-tie")

If a domain's whole surface is used as the secondary side for *multiple,
separate* tied contacts (one per target it touches -- a reasonable and
common pattern when a fan/patch-shaped domain borders more than one solid),
any node genuinely at or near the anatomical seam where those two targets
meet can be within `max_distance` of *both* primary surfaces simultaneously.
It then gets pulled toward two independent targets by two independent
contacts, and as those targets deform even slightly differently from each
other, the conflict grows over time rather than resolving. Symptom: a
domain that has multiple separate tied contacts to different targets shows
instability concentrated right at the geometric boundary between those
targets' extents, even though each contact tested alone (against a domain
that only touches one target) is stable.

Fix: track, per duplicated/tie node, which *single* target it was actually
assigned to during setup (from whatever original classification decided
"this node ties to solid X"), then when building each target's own
secondary-surface facet list, exclude any facet that contains a node
assigned to a *different* target. This preserves the "whole domain as
secondary, let distance-based engagement decide" design (gotcha 5's
workaround doesn't need to change) while preventing the specific
double-assignment case. This is not a rare edge case in practice -- on one
real model, this exclusion removed 39-69 double-tied facets per affected
domain, not just one or two isolated nodes.

### 14. Facet winding consistency needs edge-adjacency propagation, not a per-element check

If the extracted surface came from a solid element's boundary faces (hex8 or
similar), check the face-extraction table itself first (`geometry-cross-referencing.md`
Rule 9) -- a wrong or inconsistent node-index convention there produces
backwards-wound facets before winding-consistency propagation even enters
the picture, and no amount of the flood-fill below fixes a face table that's
wrong in the first place, only facets that are individually fine but
inconsistent *with each other*.

A per-facet heuristic ("does this face's normal point away from its own
parent element's centroid") corrects each facet *individually* but does not
guarantee *global* consistency across a whole extracted surface -- two
neighboring facets, each individually "outward" by that test, can still
disagree with each other if their parent elements happen to have opposite
internal vertex handedness (plausible in an unstructured/anatomical mesh
assembled from separately-generated regions). The actual property a
consistently-oriented manifold surface requires is: every shared edge is
traversed in *opposite* directions by its two owning facets. Verify this
directly (count each undirected edge's two directed traversals; flag any
edge where both facets traverse it the same way) rather than trusting a
per-element heuristic. To *fix* it, not just detect it: flood-fill from one
arbitrarily-oriented seed facet across the shared-edge adjacency graph,
flipping any neighbor whose shared edge runs the same direction as the
current facet's -- this guarantees full consistency across one connected
surface component regardless of how many individual elements have
"backwards" internal ordering.

Caveat: even after this fix, a stricter external checker (e.g. FEBioStudio's
own preprocessor) may still flag isolated facets as "incorrect winding" --
in one real case this correlated with hex elements that expose *more than
one* boundary face to the surface (a corner/edge element of the solid), and
was confirmed cosmetic (no effect on the actual solve) specifically for
`tied-node-on-facet` per gotcha 11's winding-insensitivity, but the root
cause of that residual discrepancy was not fully identified -- don't assume
it's automatically cosmetic for a different contact type without the same
kind of empirical A/B verification.

### 15. Fixing all of a shell domain's translational DOFs does not fully immobilize it

`<bc type="zero displacement">` with `<x_dof>`/`<y_dof>`/`<z_dof>` only
constrains a shell's *midsurface* translation. A shell element formulation
carries an additional through-thickness/director degree of freedom that
this BC does not cover, so a shell domain with every node "fully fixed"
this way can still deform through its thickness and report a negative
jacobian -- confirmed by fixing every node of a specific shell domain,
verifying independently that every element's node set was indeed 100%
inside the fixed `NodeSet`, and still seeing negative-jacobian reports on
that exact domain later in the same run. If the intent is to make a domain
genuinely inert (no possibility of it ever reporting a jacobian error, not
just "very stiff"), fixing its nodes is not sufficient -- delete its
`<Elements>` block entirely, **and** the corresponding `<SolidDomain>`/
`<ShellDomain>` entry in `<MeshDomains>` (leaving a domain reference with no
matching `<Elements>` block is itself an error). Keep the domain's *nodes*
declared (with the same fixed BC, if something else still needs to anchor
to them via shared node IDs) -- only the element/domain declarations need to
go. This is a different mechanism from gotcha 8's "phantom node with zero
material stiffness" near-singular-DOF problem (that one is about a node
having *no* stiffness resisting it at all; this one is about a node having
real stiffness in-plane but an *uncovered* DOF out-of-plane) -- both are
"some DOF isn't actually constrained the way it looks," but they need
different fixes, so don't assume solving one solves the other.

### 16. Converting an Ogden material's parameters between Abaqus and FEBio conventions

Both use the same underlying Ogden functional form but different symbol/
scaling conventions for the same physical model. Abaqus (and the classic
Ogden 1972 formulation): `W = Sigma (2mu_i/alpha_i^2)(lambda1^alpha_i + lambda2^alpha_i + lambda3^alpha_i - 3)`. FEBio:
`W = Sigma (c_i/m_i^2)(lambda1^m_i + lambda2^m_i + lambda3^m_i - 3)`. Matching term-by-term with
`m_i = alpha_i` gives the direct conversion: **`c_i = 2mu_i`**, **`m_i = alpha_i`**. This
also matches what FEBioStudio's GUI shows for an Ogden material's property
panel (`c1/m1`, `c2/m2`, ... up to 6 terms, plus `k`/`pressure_model`/
`density`) -- confirmed against that panel's field names directly, though the
underlying raw `.feb` XML tag/attribute names for an Ogden `<material>` block
have not been separately verified by running one through `febio4.exe` in a
minimal test file (gotcha 6) -- do that before trusting a hand-written Ogden
material block on a real model, the same as any other unfamiliar construct.

FEBio's Ogden also needs an explicit bulk modulus `k` for the volumetric
(near-incompressibility) part, which Abaqus's `poisson=` on the
`*Hyperelastic` line specifies differently; derive it from the *initial*
shear modulus and the original Poisson's ratio: `k = 2*mu0*(1+nu)/(3*(1-2nu))`.

**Correction (an earlier version of this note had the wrong `mu0` formula):**
the initial/small-strain shear modulus is simply **`mu0 = Sum_i mu_i`** -- no
factor of 2, and no dependence on `alpha_i` at all. This falls out of
differentiating the uniaxial nominal-stress formula
`sigma_nom(lambda) = (2*mu/alpha)*(lambda^(alpha-1) - lambda^(-alpha/2-1))` at
`lambda=1`: `d(sigma_nom)/d(lambda) = (2*mu/alpha)*[(alpha-1) - (-alpha/2-1)]
= (2*mu/alpha)*(3*alpha/2) = 3*mu`, independent of `alpha` -- so Young's
modulus `E = 3*mu` for a single term, and since `E = 3*G` for an
incompressible material, the term's shear-modulus contribution is exactly
`mu` (no extra factor). Summing over terms (strain energy, hence its
derivative, is additive) gives `mu0 = Sum_i mu_i` for an N-term model. A
previously-considered `mu0 = mu*alpha/2` formula is wrong in general -- it
only happens to agree with the correct one at the neo-Hookean special case
`alpha=2` (where `alpha/2=1`), which is easy to mistake for a working general
formula if that's the only case checked. This matters directly when
re-fitting a Marlow (`abaqus-inp-format.md`) or other data-driven Abaqus
material to FEBio's closed-form Ogden (`scripts/fit_ogden.pl` in this skill
does this fit), since the fit produces `(mu_i, alpha_i)` pairs in the Abaqus
convention that then need both conversions -- `c_i`/`m_i` for the deviatoric
terms, `mu0` for `k` -- before use.

### 17. `tied-node-on-facet` closes an initial gap -- it does not preserve it

Abaqus `*Tie, adjust=no` keeps whatever offset exists between the tied points.
FEBio's `tied-node-on-facet` does not: every secondary node that engages (is
within `max_distance`) is pulled onto the primary surface. Verified in a
20-line model: a free shell square 0.5 above a fixed one, tied, **no load at
all** -- the first iteration moves every node the full 0.5, the second
iteration diverges (residual 6e12, negative jacobian). With zero gap the same
file runs cleanly. What changes the outcome:

- `<node_reloc>1</node_reloc>` relocates the secondary nodes onto the primary
  surface at initialization, so the tie starts at zero gap: clean run. It moves
  geometry by up to the gap, so only use it where that shift is negligible.
- `special=0` does not help. A very low `penalty` "converges" to garbage (the
  free part drifts and rotates several mm).

This explains two otherwise-baffling symptoms: a huge `INITIAL` displacement
norm the moment a tie engages, and the knife-edge `max_distance` behaviour
(0.22 works, 0.25 fails): a larger cutoff engages nodes that are further off
the surface, and each one gets snapped. Measure the true point-to-facet gaps
(`scripts/nearest_facet_distance.pl`) before choosing `max_distance`, and use
gotcha 18 when the source intent is a node-to-node tie.

**Its `penalty` is a per-node spring in N/mm (tie force = penalty x gap), with no
auto-scaling.** Check it against the forces the tie must carry. A converted
reference model had `penalty` 0.0005 on the ties holding two lofts to the vaginal
wall. In the runs the tied points drifted a median 2.2-2.6 mm from the tissue,
which itself moved 3-6 mm, so those lofts barely held while their mirror twins
(shared nodes) held fully. A mini test pulled a loft-like strip (E 21, 0.49 mm)
0.5 mm off a soft block. The gap was 0.50 mm at penalty 0.0005 (fully detached),
0.13-0.18 at 1, 0.02 at 10, 0.003 at 100 and 0.0003 at 1000, at ~0.3 N per node.
Measure tie gaps in a run's results (tied node vs the tissue point it started on)
the way you measure any other load path. Tighter isn't always better when the tied
surfaces fold. In a three-shell folding test (shells E t = 100), penalty 1 (~1 %
of E t) ran cleanly, with the tied edges giving 1.6-2.3 mm while parts moved
35-61 mm. Penalty 100 (= E t) failed at 37 % of load, and 0.0005 tore the strip off.
That warning did not carry over to the full model, whose tied lofts fold little.
Penalty 1 and 100 (0.1x and 10x the loft's E t of 10.3 N/mm) converged alike (96 vs
99 failed attempts to full load). The tied tissue stayed within median 0.01-0.02 /
max 0.29 mm at 1, and within 0.003 mm at 100. Pick the penalty for the gap you can
accept at the expected nodal force, and test before assuming it limits convergence.

### 18. Exact node-to-node ties: `<Constraints>` linear constraints

The FEBio equivalent of an Abaqus node-to-node `*Tie` (or a tie to a
node-based master such as truss nodes) is a set of linear constraints
`u_slave - u_master = 0` per DOF. Verified syntax (FEBio 4.13, spec 4.0), as its
own top-level section after `<Contact>`:

```xml
<Constraints>
  <constraint name="LA_truss_ties" type="linear constraint">
    <tol>0.01</tol>          <!-- augmentation tolerance (relative) -->
    <penalty>10</penalty>    <!-- stiff relative to the parts being tied -->
    <maxaug>10</maxaug>
    <linear_constraint>
      <node id="22488" bc="x">1</node>
      <node id="23677" bc="x">-1</node>
    </linear_constraint>
    <!-- ...one <linear_constraint> per node pair per dof (x, y, z)... -->
  </constraint>
</Constraints>
```

- **Rejected placements** (`unrecognized tag`): `<linear_constraint>` inside
  `<Boundary>` (the FEBio 2.5 style), and `<bc type="linear constraint">`.
- A constraint whose partner node is fixed by a BC works: the free node is
  held (z = 9e-9 in the test). No special-casing of pinned nodes is needed.
- It is augmented Lagrangian: the log prints `augmentation #` lines per step.
  2-3 augmentations per step is healthy. Hitting `maxaug` every step with a
  multiplier that keeps growing by a similar fraction each time means the
  constraints are fighting a much stiffer load path (in the reference session:
  a part also rigidly merged into stiff shells). Treat that as a modelling
  conflict, not a tolerance to loosen.
- No facets, no phantom nodes, no contact search, no gap snapping. On a full
  model with 199 tied LA nodes the constraint violation stayed around 1e-6 mm
  and the tied truss chains moved plausibly (symmetric left/right).

### 19. Thick shells: averaged nodal normals can pinch elements before any load

FEBio shells are solid-like. Each element integrates through its thickness
along nodal directors that, by default, average the normals of the shell
elements around each node (`shell_normal_nodal` = 1). Where the mesh folds
more tightly than about half the shell thickness, the offset faces cross over
and the element starts out (nearly) inverted. The model still starts. It then
fails at a small load fraction with negative jacobians in the same few
elements, whatever you change about materials, contacts or solver settings. Abaqus
conventional shells (S4R/S3R) do not form this through-thickness volume, so a
thick shell imported from Abaqus can carry the defect unnoticed.

- **Detect:** `scripts/shell_pinch_check.py model.feb --domains <thick shells>`.
  Reference case: a 4 mm levator-ani shell with element edges of about 1.5-4 mm
  had 6 elements face-inverted at rest and 12 more pinched below 25 %. These
  were exactly the elements that inverted in every failing run. **Get the
  offset geometry right or this undercounts badly:** the script places the
  offset faces at 0/-t from the nodes by default (FEBio 4's actual convention
  -- nodes are the shell's TOP face, gotcha 21 below), auto-detecting
  `<shell_formulation>` for the legacy +-t/2 mid-surface case. An earlier
  version of this script always used +-t/2, which is silently wrong for a
  normal FEBio 4 model: on a second thick shell (0.49 mm, tighter folds than
  the 4 mm case above) it reported only 1 Gauss-point-inverted element where
  the correct geometry found **86**, in the same domain that kept failing at
  the same wall no matter what else changed -- don't trust an old "looks
  fine" pinch-check result without confirming which convention it used.
- **Fix:** `<shell_normal_nodal>0</shell_normal_nodal>` inside the
  `<ShellDomain>` (element normals; thickness unchanged). This took that model
  from stalling at t~0.2 to running retry-free to t~0.56, where it met an
  unrelated wall shared with the model without the part. A *different* thick
  shell fixed the same way (0.49 mm, 108 face-inverted / 86 Gauss-point-inverted
  under the correct geometry) took its model from a hard wall at t~0.86 to
  t~0.97 (99.7 % of the intended load).
- **Did not help** (tested one at a time): the Abaqus-consistent bulk modulus,
  `type="three-field-shell"`, a thinner shell with membrane-equivalent
  stiffness (worse), implicit dynamics, and `type="elastic-shell-old"` (the
  pre-2.6 mid-surface formulation -- it has no `<shell_thickness>` or
  `<shell_normal_nodal>` tag at all on `<ShellDomain>`, so thickness has to
  come from `<MeshData><ElementData type="shell thickness">` instead, and
  because it has no nodal-normal-averaging option to turn off, it re-pinches
  the same way; not a viable escape hatch for this defect).
- Other `<ShellDomain type="...">` values in FEBio 4.13 (from `febiomech.dll`,
  gotcha 7): `elastic-shell`, `elastic-shell-ans`, `elastic-shell-eas`,
  `elastic-shell-old`, `three-field-shell`, `rigid-shell`. `three-field-shell`
  needs an uncoupled material ("Model initialization failed" with
  neo-Hookean, runs with the uncoupled `Ogden`). `elastic-shell-eas` failed a
  small Ogden test that the default shell passed.

### 20. Plot output: per-iteration plotting, file size, and "converged" flags

- `<plot_level>PLOT_MINOR_ITRS</plot_level>` writes every Newton iteration,
  including the failed ones after the last converged step. That is exactly
  what you need to see where a blow-up starts (`scripts/feb_postmortem.py`).
  It is large: a 55-step run with ~400 failed iterations produced a 2.6 GB
  `.xplt`. `<plot_range>` did not limit it in practice, so check free disk
  space first.
- FEBioStudio "debug" runs flag converged states with status 0 and iteration
  states with 2. A plain `febio4 -i` run with `PLOT_MINOR_ITRS` flags *every*
  state 2, so identify converged states from the log's
  `converged at time` lines (`xplt_reader.converged_state_indices`).

### 21. A shell contact/tie needs an `offset` when the Abaqus source shell was meshed mid-surface

Abaqus S4R/S3R shells are meshed at the geometric mid-surface: the physical
shell occupies `[nodes - t/2, nodes + t/2]` along the normal, and an
`*Elset`-based `SPOS`/`SNEG` surface used in a `*Contact Pair` sits exactly
`t/2` off the nodes. FEBio 4's default shell geometry is different (gotcha
19): the nodes themselves ARE the top face, and the shell occupies
`[nodes - t, nodes]`. A FEBio `<contact>`/tie built straight from the shell's
own `<Surface>` (built from its element facets, i.e. the node positions) acts
at the nodes, not at the mid-surface -- so a contact converted node-for-node
from an Abaqus `SPOS`/`SNEG` pair engages `t/2` too early or too late,
compared to the source.

Verified in a 20-line model (a hex block over a 4-unit-thick shell plate,
pushed down): with `<offset>0</offset>` on a `sliding-elastic` contact, the
force onset exactly tracked the node surface (no allowance for the shell's
own half-thickness); with `<offset>2</offset>` (half the 4-unit thickness),
onset shifted by exactly that amount, matching a mid-surface-meshed Abaqus
source. `offset` works this way in both `two_pass="0"` and `"1"` mode.

**Fix:** set the contact/tie's `<offset>` to the Abaqus shell's half-thickness
(the side being approached determines the sign convention Abaqus itself would
use for SPOS vs SNEG; get this by comparing where the counterpart surface
actually sits relative to the shell's mid-surface, not by guessing). Reference
case: an Abaqus frictionless contact between a solid's back face and a full
levator-ani shell surface (`SPOS`), converted as FEBio `sliding-elastic` with
`<offset>2</offset>` on a 4 mm shell, reproduced the source geometry exactly
(confirmed by measuring the true point-to-facet gap on both sides, coordinate
by coordinate) -- even though in that particular case the gap never closed
enough during the run for the contact to actually engage.

### 22. `<discrete dmat="k">` resolves the k-th `<discrete_material>` by LIST POSITION, not by its `id` attribute

Deleting a `<discrete_material>` block (e.g. removing some stand-in ground
springs that turned out not to be in the Abaqus source) and leaving the
remaining blocks' `id` attributes as they were produces
`tag "discrete" (line N) : invalid value for attribute "dmat"` on every
`<discrete dmat="...">` binding that comes after the gap, even though the
referenced `id` still exists elsewhere in the file and looks perfectly valid
read in isolation. FEBio does not look up `dmat` against the `id` attribute
at all -- it indexes into the `<discrete_material>` elements in file order,
1-based. Deleting the 2nd of 5 materials means what used to be `dmat="3"` is
now the *2nd* material in the list, so any `<discrete>` that still says
`dmat="3"` now resolves to the wrong material (or, past the end of the list,
errors outright).

**Fix:** after adding or removing any `<discrete_material>`, walk every
remaining one in file order, renumber its `id` to `1..N`, and rewrite every
`<discrete dmat="...">` binding to match the new number for that same
material -- don't just leave old `id` values in place and assume they still
resolve correctly. This mirrors gotcha 1's node-ID lesson: FEBio's on-disk
`id`/index conventions are stricter and more positional than they look.

### 23. A Marlow material refit must be I1-only (use `Yeoh`), and a constant `k` from mu0 softens it at large strain

Abaqus `*Hyperelastic, marlow` builds its deviatoric energy from one test
(usually `*Uniaxial Test Data`) as a function of **I1 only**. That fixes its
biaxial and planar response too. A 1-term Ogden fitted to the same uniaxial
points is a different model once you leave uniaxial tension. For a soft
muscle curve (0.0058 MPa at 10 %, 0.2963 MPa at 100 %), an Ogden fit with
c1=0.0427, m1=6.49 was within about 10 % of the data in uniaxial tension but
**1.4x softer in equibiaxial stretch at 10 %, 2.1x at 20 %, 2.6x at 30 % and
2.9x at 40 %**. Planar stretch was 1.2-1.35x softer. For a large exponent,
Ogden's biaxial response is barely stiffer than its uniaxial response. For an
I1 law, equibiaxial stretch l has the I1 of a much larger uniaxial stretch.
A pressure-loaded sheet is in biaxial tension, so this matters exactly where
the refit is used.

- **Use FEBio's uncoupled `Yeoh`** (`<c1>`..`<c6>`, `<k>`, confirmed in 4.13;
  also listed in `febiomech.dll`: `polynomial`, `uncoupled isotropic
  Lee-Sacks`, generic `hyperelastic`). Yeoh is I1-only, and its uniaxial
  nominal stress `2(l - l^-2) Sum i c_i (I1-3)^(i-1)` is linear in `c_i`.
  The fit is a unique linear least-squares solve with no multi-start. In the
  case above, 2 terms (c1=0.0100393, c2=0.0172014, both positive, hence
  stable) matched Marlow within about 10 % in uniaxial, planar and
  equibiaxial stretch. `scripts/fit_yeoh.py` does the fit and the
  three-state comparison. `mu0 = 2 c1`.
- **`k` from mu0 only matches nu at zero strain.** `k = 2 mu0 (1+nu)/(3(1-2nu))`
  (gotcha 16) is constant in FEBio. As the deviatoric response stiffens about
  10x, the material becomes relatively compressible: J = 1.06-1.23 at 50 %
  strain in single-hex tests. That costs 20-45 % of the stress. Abaqus Marlow
  with `poisson=0.47` keeps nu constant (uniaxial J = l^(1-2nu), 1.025 at
  50 %). Matching that at 20/30/40/50 % uniaxial strain needed k =
  0.52/0.78/0.99/1.33 MPa versus 0.33 from mu0. `fit_yeoh.py` prints these
  values. Pick k for the working strain range and treat it as its own
  single-variable change.
- **Expect a modest effect on a pressurised membrane's deflection.** Bulge
  scales roughly with the cube root of membrane stiffness. In a pinned
  60x60x4 mm panel test, Ogden -> Yeoh -> Yeoh with k x10 moved the centre
  deflection 17.6 -> 16.4 -> 15.9 mm at 20 % of the pressure and 42.6 -> 36.1
  -> 29.8 mm at full pressure. It is the faithful conversion, but it will not
  remove a large bulge that the source material itself produces.
- Verification pattern: single hex8 (10x10x4) with symmetry BCs for uniaxial,
  planar and equibiaxial stretch. The nominal stress is `sx*J/l` from
  `<logfile><element_data data="sx;sy;sz;J">` (`J` is a valid log variable).
  The FEBio Ogden hex matched the analytic formula, which confirms gotcha 16's
  `c_i = 2 mu_i`.

### 24. Single-element shell tests: pull the back face too, or the shell looks 2-3x too soft

A FEBio 4 shell node carries the top-face displacement (x/y/z) *and* the
back-face displacement (`sx`/`sy`/`sz`). Prescribing only x/y/z on a one-quad
shell (10 mm, 4 mm thick) leaves the back face free to lag. The reported
stress was **2-2.7x lower** than the same material in a hex8, and nearly the
same in uniaxial, planar and equibiaxial stretch. Prescribing the back face
in-plane as well made the shell match the hex8 to 4 digits. Leave `sz` free so
the thickness can change. Syntax (FEBio 4.13, verified):

```xml
<bc name="x1_sx" node_set="x1" type="prescribed shell displacement"><dof>sx</dof><value lc="1">5</value><relative>0</relative></bc>
<bc name="x0_sx" node_set="x0" type="zero shell displacement"><sx_dof>1</sx_dof><sy_dof>0</sy_dof><sz_dof>0</sz_dof></bc>
```

This is a test-setup trap, not an element defect. A pressurised pinned
60x60x4 mm panel (12x12 quads) deflected the same as a 4-layer hex8 of the
same body with the same support: 17.6 vs 17.6 mm at 20 % pressure, and 42.6
vs 43.3 mm at full pressure. Support matters only at large load. Pinning just
the node (top) face instead of the whole thickness added about 13 % at full
pressure and about 1 % at 20 %. That is a real 3D effect, and the hex shows it
too. An Abaqus S4R pinned at mid-surface nodes transmits the support through
the whole thickness.

Also verified: a two-point `SMOOTH STEP` load curve starting at `t0 > 0`
(`<pt>0.5,0</pt><pt>1,1</pt>`, `extend CONSTANT`) holds 0 until t0 and then
ramps. This is the way to put one load on a later ramp.

Also verified (4.13): **truss elements, e.g. to give massless spring chains
their Abaqus truss mass in a DYNAMIC run.** Write them as `<Elements type="line2">` with
`<BeamDomain name=".." mat=".." type="linear-truss"><cross_sectional_area>1</cross_sectional_area></BeamDomain>`
and a `<material type="linear truss">` with `<density>`, `<E>` and `<v>`. The solver also accepts the older form,
`<Elements type="truss2">` under a `<SolidDomain>` with the same `cross_sectional_area`, and gives identical results
(same displacement to 12 digits in a two-node test). But **FEBio Studio 3.2 cannot open a file that uses
`truss2`**: "tag "Elements" : invalid value for attribute "type"". Its reader knows `line2`/`line3`,
`BeamDomain` and `linear-truss`. So use the `line2` form, or convert old files with a text-level rewrite of just
those lines (`claude_diag`-style `feb_truss2_to_line2.py`). Without the area FEBio stops with "Cross sectional area
of "X" must be positive", and `<area>` is an unrecognized tag. A two-node test with a step force gave the
consistent mass rho*A*L/3 at the free node. With `E` ~ 1e-6, such trusses add mass and almost no stiffness.

Also verified: **mass-proportional damping** is a body load,
`<Loads><body_load name="md" type="mass damping"><C lc="3">20</C></body_load></Loads>`.
`C` is in 1/s and the damping force is -C M v. A two-node test gave 10.02 mm against the analytic 10.03 mm
for C = 10. Put a load curve on `C` so the damping switches on only after the loading ramp (dynamic
relaxation; see convergence-debugging.md).

### 25. Beams (`linear-beam`): the `BeamDomain` needs a `type`, the element shear-locks, and two pins leave a mechanism

FEBio 4.13 has a geometrically exact (Simo-Reissner) elastic beam. Verified syntax:

```xml
<material id="9" name="bm" type="linear-beam"><density>1e-9</density><E>10</E><G>4</G>
  <A>0.01</A><A1>5</A1><A2>5</A2><I1>4</I1><I2>4</I2></material>
<Elements type="line2" name="B"> ... </Elements>
<MeshDomains> <BeamDomain name="B" mat="bm" type="linear-beam"/> </MeshDomains>
<bc name="clamp" node_set="root" type="zero rotation"><u_dof>1</u_dof><v_dof>1</v_dof><w_dof>1</w_dof></bc>
```

- `A` is the axial area, `A1`/`A2` the shear areas and `I1`/`I2` the bending moments of inertia. Torsion
  uses J = I1 + I2 (an L-shaped cantilever matched to 0.2 %). A and I are independent, which a shell tube
  can't offer (EI/EA = R^2/2).
- Without `type="linear-beam"` the parser says `tag "BeamDomain" : unrecognized tag`, and `SolidDomain` is
  rejected for `line2`. The parameter names come from the SDK header
  (`FEBioStudio/sdk/include/FEBioMech/FEElasticBeamMaterial.h`) and the domain syntax from the shipped
  manual (`FEBioStudio/doc/FEBio_User_Manual.pdf`, section "BeamDomain"; `pdftotext` from the MinGW
  toolchain extracts it). Check both before guessing.
- **The 2-node element is fully integrated, so it shear-locks.** Each element bends with
  EI + G*A_s*L^2/12. A 10-element cantilever gave 0.2530 mm where plain Timoshenko theory gives 0.3358 and the
  locked formula gives 0.2525. Consequence: a two-element zig-zag of a beam chain is resisted by G*A_s,
  not EI (G*A_s is the short-wave buckling load). Tiny shear areas make the beam useless against it.
- A curved beam pinned only at its two ends can swing rigidly about the chord. The first iteration then
  has a displacement norm around 1e17. A straight one with its twist unloaded runs. Restrain that mode
  (end rotations, or ties to the surroundings).
- The beam's rotations appear as extra equations (+3 per beam node, minus those fixed). Compare the log's
  "Nr of equations" with and without the beam to confirm they are active.

**It did not work as a stabiliser in the reference model.** A beam in parallel with the 18-spring
posterior-arcus chains, whose nodes are tied to the LA by linear constraints, was supposed to add
bending/shear stiffness with almost no axial stiffness. Every variant broke the static run earlier than
without it: stiff (G*A_s 20 N, EI 40) with clamped ends failed at t = 0.018, 10x softer at 0.065, stiff with
pinned ends at 0.13, and a non-symmetric global matrix at exactly the same 0.13. The diverging attempts ran
away at an LA-perineal body connector region (negative jacobians in the PeB and LA), not at the arcus.
The cause is not yet known. Test the combination (beam + linear-constraint ties + tension-only connectors)
in a mini model before relying on it.
