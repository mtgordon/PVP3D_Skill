# Repairing degenerate/sliver mesh regions with mmg (mmgs)

Sometimes a specific shell domain has genuinely broken geometry -- a sliver
triangle with nearly-collinear nodes, baked into the mesh itself rather than
caused by deformation -- and no amount of contact/solver tuning will fix a
negative jacobian coming from an element that starts out degenerate. This is
a companion technique to `convergence-debugging.md`: use it once the
diagnostic loop there (specifically the interactive-console per-element
report) has pointed at a *specific, consistently-failing element ID* that
turns out to have a near-zero area or a near-180 deg/near-0 deg interior angle
even in the undeformed reference configuration (Rule 4 in
`geometry-cross-referencing.md` gives the area/aspect/angle formulas to
confirm this). `mmgs` (the surface-remeshing tool from the mmg project) can
clean up that region's interior mesh quality while leaving its boundary
exactly where it was.

## Getting mmgs without admin rights

`mmgs` can be built from source with a fully portable toolchain -- no
installer, no admin privileges needed:

- A portable MinGW-w64 GCC distribution (e.g. WinLibs) for the compiler.
- A portable CMake (the `.zip` distribution, not the installer).
- Build mmg from its source release with these pointed at via
  `-DCMAKE_C_COMPILER=...`/standard CMake generator flags; the resulting
  `mmgs.exe` (surface remesher -- there's also `mmg3d` for volume meshes and
  `mmg2d`, not needed for a shell/surface repair) is a standalone
  executable with no further runtime dependencies.

## The core technique: preserve the boundary, remesh only the interior

**The critical design choice, learned from a real failure**: do not let
mmgs freely remesh a patch's entire extent including its boundary, even
though that produces the visually "cleanest" output. A patch's boundary
nodes are usually where it anatomically/structurally attaches to
neighboring geometry (shared nodes with an adjacent part, or the specific
points a tied contact targets) -- regenerating them lets the remesher
reconstruct the attachment path using *only the patch's own, locally
degenerate* input geometry, disconnected from the true position of
whatever it's supposed to connect to. The result can be numerically
stable (good element quality throughout) while being **structurally/
anatomically wrong** -- connecting to the wrong location entirely. This
class of bug is easy to miss from log files and solver behavior alone
(nothing about it causes a convergence failure) and was only caught by a
human visually inspecting the deformed output shape -- a concrete instance
of `convergence-debugging.md`'s "when to stop guessing and ask for human
input" advice, except here the trigger wasn't a stuck solver, it was a
result that looked fine numerically but wrong anatomically.

The fix: mark **every** original boundary node -- both nodes shared with
other domains (fan-to-fan/fan-to-tube style) and nodes that are the target
of some other connection (fan-to-solid tie style) -- as a `RequiredVertices`
entry in mmgs's input mesh. mmgs then only touches genuinely interior
nodes, guaranteeing every attachment path stays exactly where the source
geometry defines it, while still fixing the interior sliver triangles that
were the actual problem.

## The `.mesh` format (MEDIT ASCII), minimal working subset

```
MeshVersionFormatted 2
Dimension 3
Vertices
<N>
<x> <y> <z> <ref>          <!-- one line per vertex, ref is usually 0 -->
...
Triangles
<M>
<n1> <n2> <n3> <ref>        <!-- 1-based LOCAL vertex indices, not your
                                  original node IDs -- build a local
                                  numbering when writing this file, and
                                  keep the local->original ID mapping
                                  around to translate results back -->
...
RequiredVertices
<K>
<local-vertex-index>        <!-- one per line, K of them -- these are
                                  the ones mmgs must NOT move or remove -->
...
End
```

Run: `mmgs.exe -in patch.mesh -out patch_remeshed.mesh -hmax <h> -hmin <h>
-hausd <h> -hgrad 1.3 -nr -v 0`. The output `.mesh` has the same format;
read its `Vertices`/`Triangles` blocks back and re-map local indices to
whatever ID scheme the target format needs (see `abaqus-inp-format.md`/
`febio-xml-format.md` for placing new nodes correctly -- remember gotcha 1's
node-ID-monotonicity rule if the target is a `.feb`).

## Choosing `-hmax`/`-hmin`/`-hausd` without guessing

Derive these from the patch's own existing geometry rather than picking
arbitrary numbers: compute every edge length across the patch's current
triangles, take the median as a characteristic length `L`, then a
reasonable starting point is `hmax ~ 1.5*L`, `hmin ~ 0.15*L`,
`hausd ~ 0.05*L` (the Hausdorff-distance parameter controlling how far the
new surface may deviate from the input -- small relative to `L` keeps the
remeshed surface close to the original shape), with `-hgrad 1.3` limiting
how fast element size can change from one region to the next. This scales
automatically across differently-sized patches in the same model instead of
needing per-patch manual tuning, and was validated across multiple
differently-shaped/sized regions in the same project without further
adjustment.

## Boundary node handling after remeshing: match on coordinate, not order

mmgs's output vertex order does not correspond to the input order, even for
`RequiredVertices` (which keep their *position* exactly, but not
necessarily their *index*). After remeshing, match each output vertex back
to its role by exact coordinate (Rule 3 in `geometry-cross-referencing.md`):
- A vertex whose coordinates exactly match an original *shared* boundary
  node -> reuse that original ID verbatim (it's still the same shared node).
- A vertex whose coordinates exactly match an original node that's the
  target of a *tie* elsewhere -> this is where the tie-duplication pattern
  from `febio-xml-format.md` gotchas 11-13 applies; treat it the same way
  the pre-remesh version did (duplicate + tie, or shared, depending on what
  that specific boundary node is for) rather than assuming the remesh
  changes the connection strategy.
- Everything else is a genuinely new interior vertex -- assign it a fresh ID
  per the target format's numbering rules.
