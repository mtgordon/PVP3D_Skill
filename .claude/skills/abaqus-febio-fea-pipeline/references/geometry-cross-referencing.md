# Working with multi-megabyte FE model files

Both `.inp` and `.feb` files in this domain are routinely tens of thousands of
lines and several megabytes. Every technique below exists because the naive
approach (`Read` the whole file, or a full nested-loop comparison) either
blows the context window or is too slow to iterate with.

## Rule 1: never load the whole file into context

Use `grep`/`sed`/`awk` to pull out only the section you need. A file-reading
tool's line-range `offset`/`limit` is fine for *looking at* a known location,
but *finding* that location, and *comparing* two files, should go through
shell tools. Useful patterns (bash, but the same ideas port to PowerShell):

```bash
# Locate a named block and print just its extent (works for FEBio XML blocks
# and Abaqus keyword blocks alike, since both use a clear start/end marker)
awk '/<Surface name="X">/{flag=1} flag{print} /<\/Surface>/{if(flag)exit}' file

# Get just the line numbers, cheaply, before deciding whether to read/edit
awk '/<Surface name="X">/{s=NR} /<\/Surface>/{if(s){print s","NR; exit}}' file
```

## Rule 2: build an ID->coordinate lookup once, reuse it for every query

Comparing node positions between an Abaqus part (part-local IDs) and its
merged/renumbered FEBio equivalent, or checking distances between two named
node groups, is fundamentally a lookup problem. Don't re-grep the whole file
per node -- that's O(n*m) against a huge file and gets slow fast. Build one
flat `id|x,y,z` table with a single pass (`scripts/build_node_lookup.sh` in
this skill does this), then do all subsequent lookups in-memory via awk
associative arrays (`scripts/nearest_point_distances.sh` does the common
nearest-point query):

```bash
# One-time extraction (works for FEBio's <node id="N">x,y,z</node> lines)
grep -oE '<node id="[0-9]+">[^<]+</node>' model.feb |
  sed -E 's/<node id="([0-9]+)">([^<]+)<\/node>/\1|\2/' > all_nodes.txt

# Nearest-point / minimum-distance query against a reference point set,
# reusing the lookup table built above
awk -F',' -v reffile="reference_points.txt" '
BEGIN{
  n=0
  while ((getline line < reffile) > 0) { split(line,a,","); n++; px[n]=a[2]; py[n]=a[3]; pz[n]=a[4] }
}
{
  x0=$2; y0=$3; z0=$4; mind=1e9
  for(i=1;i<=n;i++){ dx=px[i]-x0; dy=py[i]-y0; dz=pz[i]-z0; d=sqrt(dx*dx+dy*dy+dz*dz); if(d<mind) mind=d }
  print $1, mind
}' target_points.txt
```

A very common mistake here: extracting a `<...>` attribute's numeric ID with a
naive regex also matches unrelated `id="N"` attributes nearby (e.g. an
element's local `id=` inside a `<Surface>` block, which restarts at 1 and has
nothing to do with the node IDs in its connectivity list). Always anchor the
extraction to the specific tag (`<node id=`, not just `id=`), and when
extracting connectivity, strip the tag around the comma-separated values
before splitting, not the raw line.

## Rule 3: match geometry by coordinate, never by ID, across a conversion

Abaqus part-local IDs and FEBio's flat renumbered IDs have no relationship.
The only reliable cross-reference is exact (or near-exact, allowing for tiny
export round-tripping differences) coordinate match. Round coordinates to a
sane number of decimals before comparing if two independently-exported files
might differ in the last couple of significant figures; but note that nodes
genuinely coincident *within the same file* (e.g. because two originally
separate parts were welded together during meshing) usually agree to full
float precision, since they came from the same source geometry -- a
loose-tolerance match across *different* files but a tight-tolerance match
*within* one file are different, both-useful checks.

```bash
# Build reverse lookup: coordinate -> node id(s), for the "other" file/region,
# then check which of "my" nodes have an exact match
awk -F'|' '
BEGIN{ while((getline line < "other_nodes.txt")>0){split(line,a,"|"); seen[a[2]]=seen[a[2]]" "a[1]} }
{ if($2 in seen) print $1, seen[$2] }
' my_nodes.txt
```

Two distinct outcomes both matter and mean different things:

- **A single shared node ID used by both regions** (i.e. the "same" node
  literally appears once, referenced by elements from both named parts) --
  this means the mesh was welded/merged at that point: displacement there is
  automatically shared, no contact/tie is needed or possible.
- **Two distinct node IDs at the identical coordinate** -- this means the
  regions are geometrically coincident but *not* mechanically connected
  unless something (a tie, a shared boundary condition, a discrete spring)
  explicitly links them. Finding this and mistaking it for a real connection
  is an easy way to think two parts are joined when they are not (or vice
  versa) -- always check which case you're in before drawing conclusions
  about model connectivity.

When doing this at scale (e.g. "does *any* node of part A coincide with *any*
node of the rest of the mesh?"), also filter matches against which nodes are
actually referenced by real elements -- a converted/re-exported mesh often
leaves behind orphaned "isolated vertex" leftovers (see `febio-xml-format.md`
gotcha 9) that will spuriously show up as coordinate matches to something
that is, mechanically, dead data.

The same trap applies to **reading results** at a coordinate-matched node. A
nearest-node lookup returns the lowest ID at that point, which can be an orphan
duplicate that never moves. In the reference model the Abaqus arcus coordinates
matched both an orphan node (ID ~1670) and the live chain node (ID ~23700). A
connector-elongation audit read the orphan and reported a 2.9 mm median instead
of 17-18 mm. When several nodes share a coordinate, pick the one used by an
element (or a spring), and cross-check one number against an independent run.

## Rule 4: checking element/mesh quality without a GUI

Ruling out "a badly-shaped element is causing the solver to choke" as a
hypothesis doesn't require opening a mesh viewer. For a planar-ish
quad/tri facet, compute: area (half the cross product of the diagonals, or
of two edge vectors for a tri), edge-length min/max ratio (aspect ratio),
and interior angles (via `atan2` of the cross/dot product of adjacent edge
vectors). A one-pass awk script over an element connectivity list plus the
node lookup table from Rule 2 gives area/aspect/min-angle/max-angle for
every element in the region of interest in seconds. Flag anything with a
near-zero area, an aspect ratio above ~5-10, or a minimum interior angle
under ~10-15 deg as suspect -- but expect most real biomechanical/anatomical
meshes to have *some* skewed elements already, so this is a "rule out
something exceptional" check, not a "the mesh must be perfect" check.

## Rule 5: editing a huge file -- line-range replace, and the reverse-order trap

The `Edit` tool's exact-string matching is fine for short, unique snippets,
but becomes unreliable and slow for blocks spanning hundreds of lines (easy
to introduce a whitespace mismatch, and the tool has to search the whole
file for the match). For large block replacement, it's more reliable to
work with explicit line numbers via `awk`:

```bash
# Replace lines [start, end] (inclusive) with the contents of replacement.xml
awk -v s="$start" -v e="$end" -v repl="replacement.xml" '
  NR==s { while ((getline line < repl) > 0) print line; skip=1 }
  NR==e { skip=0; next }
  !skip { print }
' file.feb > file.feb.tmp && mv file.feb.tmp file.feb
```

**If replacing multiple blocks in the same file, always work from the
highest line numbers down to the lowest.** Replacing an earlier block first
shifts every line number after it, silently invalidating the line ranges you
already computed for later blocks. Either process strictly in reverse, or
re-`grep` the current line numbers fresh before every single replacement --
don't trust line numbers computed before an earlier edit in the same file.

After any such surgery, **always validate well-formedness before spending
time running the solver on it** -- a cheap, fast check that catches a huge
class of "forgot to close a tag" / "mismatched block boundaries" mistakes
before burning a solver run on them (`scripts/validate_feb_xml.ps1` in this
skill):

```powershell
try {
  [xml]$xml = Get-Content -Raw "model.feb"
  Write-Output "XML VALID"
} catch {
  Write-Output "XML ERROR: $($_.Exception.Message)"
}
```

This validates well-formedness only, not FEBio's stricter semantic rules
(see `febio-xml-format.md`) -- a file can pass this check and still fail to
parse in FEBio. Treat it as a fast pre-filter, not a full validation.

## Rule 6: the "add a facet-only surface with no material" pattern

> **Superseded -- kept for history.** A later session found this pattern is
> singular by construction (`febio-xml-format.md` gotcha 8, resolution note).
> Tie the chain's nodes to the shell's nodes with FEBio linear constraints
> instead (`febio-xml-format.md` gotcha 18), resolving the pairs from the
> Abaqus `*Tie` (`abaqus-inp-format.md`, "Resolve every `*Tie`"). No phantom
> nodes, facets or contact are needed.

When a source model ties a bare 1D chain (no natural facets) onto a shell
surface, and the target format requires real facet elements for any tied
contact surface (see `febio-xml-format.md` gotcha 5), one workable pattern
is:

1. Keep the chain's real mechanical behavior on its own true centerline
   nodes (e.g. as a chain of discrete spring/truss elements -- the actual
   structural stiffness).
2. For each centerline node, add one "phantom" partner node offset by a
   very small, fixed vector (small enough to be geometrically negligible,
   e.g. two to three orders of magnitude smaller than the model's own
   contact tolerances -- pick a vector with nonzero components along every
   axis so it's never accidentally parallel to a chain segment's own
   direction, which would produce a degenerate zero-area facet).
3. Build a ribbon of quad facets connecting consecutive
   centerline/phantom node pairs -- this is a valid `<Surface>` for tie
   purposes even though it carries no material of its own.
4. **Do not leave the phantom nodes mechanically unsupported** -- see
   `febio-xml-format.md` gotcha 8. This step was where the reference
   session ran out of time; if picking this up again, verify the
   phantom-to-centerline connector's stiffness in the minimal test-file
   harness (with a single tied contact against a stand-in shell) before
   trying it on the full model, and check whether a node shared by
   adjacent ribbon facets is being tied redundantly.

## Rule 7: a node "not shared with any other domain" is not the same as a node "on the domain's true boundary edge"

When looking for a domain's genuine attachment/free edge (e.g. to find
which of its nodes should connect to some other structure), it's tempting
to check "which of this domain's nodes are not referenced by any *other*
domain" -- but that query also matches every ordinary *interior* node, since
an interior node is by definition never shared with anything else either.
On a domain with a large interior and a small true boundary, this
over-inclusive check can return many times more candidate nodes than
actually exist on the real edge, and the excess candidates then need
threshold-based distance filtering to weed out (unreliable -- see Rule 4's
"some skewed elements are normal" caveat, the same applies to "some
false-positive proximity is normal").

The correct check is genuine mesh-boundary-edge detection: for a domain's
own triangle/quad element list, an edge (an unordered pair of adjacent
vertex IDs) that appears in only **one** element of that domain is a true
boundary edge; every node touching at least one such edge is on the real
boundary. An edge appearing in two elements is interior to the domain and
tells you nothing about its outer boundary, regardless of what else those
two elements do or don't share with other domains.

```bash
# one-pass boundary-edge detection over a domain's own element list
# (element_nodes.txt: one element per line, space-separated vertex IDs)
awk '{
  n = NF
  for (i=1; i<=n; i++) {
    a=$i; b=$(i%n+1)
    key = (a<b) ? a","b : b","a
    count[key]++
  }
}
END {
  for (k in count) if (count[k]==1) { split(k,p,","); boundary[p[1]]=1; boundary[p[2]]=1 }
  for (nid in boundary) print nid
}' element_nodes.txt
```

Run this *before* any node-domain-sharing or distance-based filtering, not
instead of it -- the two checks answer different questions ("is this
geometrically on the edge" vs. "is this already connected to something
else") and both are usually needed together to find the specific subset of
boundary nodes that are still unconnected.

## Rule 8: measuring the real gap a tied-contact uses -- point-to-facet, not point-to-point

`nearest_point_distances.sh` (Rule 2) measures the distance from a point to
the *nearest vertex* of a reference set. A proximity-based facet contact
(e.g. FEBio's `tied-node-on-facet`) actually projects onto the nearest
*point on a facet's surface* -- which can be meaningfully closer than any of
that facet's own vertices (the closest point can land on a facet's
interior or an edge between two vertices, not at a vertex at all). Using
the vertex-based distance as a stand-in for the real engagement gap can
overstate it by a wide margin -- enough, in one real case, to make a
gap that was actually a few hundredths of a unit look like several tenths.
When the exact tolerance value matters (not just "roughly enough gap"),
use `scripts/nearest_facet_distance.pl`, which computes the true
closest-point-on-triangle distance (splitting any quad facets into two
triangles first) instead of nearest-vertex. This is also the right tool for
confirming whether a `max_distance` pass/fail boundary is driven by a
specific node crossing the threshold, or by something else entirely (see
`convergence-debugging.md`'s note on chaotic contact-parameter
sensitivity) -- compare the *set* of nodes/facets within distance at each
tested threshold, not just whether the run succeeds or fails.

## Rule 9: extracting a solid element's boundary faces needs a correct, consistent node-index convention

A hex8 (8-node brick) element's six faces are each a specific 4-node subset
of its own local node ordering -- get this table wrong (or use an
inconsistent convention across faces) and every extracted "outward" facet
normal comes out backwards, which then surfaces later as a facet-winding
bug (`febio-xml-format.md` gotchas 11/14) in whatever contact/tie uses that
extracted surface. The standard convention (0-indexed into the element's own
8 nodes, matching FEBio's and most other codes' hex8 node ordering):

```
face 0: nodes 0,1,2,3   face 1: nodes 4,7,6,5
face 2: nodes 0,4,5,1   face 3: nodes 1,5,6,2
face 4: nodes 2,6,7,3   face 5: nodes 3,7,4,0
```

To find a solid domain's true *boundary* faces (as opposed to internal faces
shared between two adjacent solid elements): enumerate all six faces of
every element using this table, key each face by its sorted node-ID tuple
(so two elements sharing a face produce the same key regardless of which
element declared it or in what rotation), and keep only the faces whose key
appears exactly once across the whole domain -- those are the true exterior
boundary, in the same "appears in only one element" spirit as Rule 8's
2D edge-boundary detection, just one dimension up. Verify the result against
a known reference (e.g. the domain's overall bounding box, or a specific
face you can independently confirm should/shouldn't be on the boundary)
before trusting it for a tie/contact surface -- an off-by-one or
inconsistently-ordered face table produces a plausible-looking but wrong
surface with no error message anywhere.

## Rule 10: always check for the user's own prior attempts first

Before re-deriving anything from scratch, look for files that already look
like earlier iterations on the same problem -- differently-named variants,
sibling `.log`/`.xplt` output next to a `.feb` with a suggestive name
(`*_best*`, `*_v2*`, `*_notubeelems*`, `*_fixed*`, `*_diag*`, etc.), or old
run logs sitting in the same directory. A quick `ls -la` sorted by
modification time, or a directory listing filtered by a shared name prefix,
routinely surfaces hard-won prior discoveries (a parameter value someone
already found by trial and error, a specific contact/BC combination that
was already known to work better) that would otherwise take another full
debugging pass to rediscover independently. Running the user's own
candidate file fresh (don't just trust a stale log sitting next to it -- logs
can predate later edits to the same file) is cheap and can reset the whole
investigation onto a much better baseline.

## Rule 11: check a thick shell's through-thickness geometry, not just its mid-surface

Rule 4's area/aspect/angle checks look at the mid-surface only. A FEBio shell
also has a thickness laid along (by default) averaged nodal normals, and where
the mesh folds tighter than about half the thickness, the element's offset faces
cross over before any load is applied (`febio-xml-format.md` gotcha 19). Rule 4
passes such elements (reference case: aspect ratio 1.5-2.2, no slivers). Run:

```bash
py -3 scripts/shell_pinch_check.py model.feb --domains ThickShellA,ThickShellB
```

It reports, per element, the offset-face and Gauss-point area ratios relative to
the mid-surface (<= 0 = inverted at rest, < 0.25 = badly pinched) and lists the
worst elements. Do this whenever the same few elements of a thick shell fail in
every variant you try. In the reference session those elements had passed
every mid-surface quality check, and they were the whole story.

## Rule 12: after any mesh surgery, audit which boundary conditions still act on anything

A BC is attached through node IDs. Remeshing (mmg), duplicating nodes for a tie,
or re-exporting from FEBioStudio can move the elements onto new node IDs at the
*same coordinates* while the BC's `NodeSet` still lists the old ones. The BC then
pins orphans, which FEBio silently drops ("N isolated vertices removed"), and
nothing warns you. This is Rule 3's "two distinct IDs at one coordinate" trap
applied to boundary conditions. Audit active-versus-declared nodes per BC
across versions:

```bash
py -3 scripts/bc_audit.py original.feb remeshed.feb current.feb
```

A count that drops (e.g. 26/26 -> 0/26) is a lost support. The script also names
the live nodes sitting on the orphaned coordinates, which are the ones to add
to the set. Reference case: two supports (a ligament fan's anchor apex and a
perineal-body anchor edge) had been silently free for every run since a remesh.
Re-pinning them turned a 79-retry run into a retry-free one up to the same
final wall.
