#!/usr/bin/env bash
# For every point in targets.txt, find the minimum distance to any point in
# reference.txt. Both files use the "id|x,y,z" or "id,x,y,z" lookup format
# produced by build_node_lookup.sh (either separator works). Useful for:
#   - measuring the real gap a tied/sliding contact has to bridge, to pick
#     an appropriate proximity tolerance (references/convergence-debugging.md)
#   - checking whether a converted surface's nodes have drifted off an
#     original Abaqus part's true centerline (references/abaqus-inp-format.md)
#   - detecting merged vs. merely-coincident nodes (references/geometry-cross-referencing.md
#     Rule 3) -- pair this with an exact-match check, since "distance ~0" and
#     "distance == 0.000000" mean different things depending on float
#     precision from independent CAD exports
#
# Usage: nearest_point_distances.sh targets.txt reference.txt
# Output: "<target_id> <min_distance>" per line, unsorted (pipe through
#         `sort -k2 -n` to see closest/farthest first)

set -euo pipefail
targets="$1"
reference="$2"

awk -F'[|,]' -v reffile="$reference" '
BEGIN{
  n=0
  while ((getline line < reffile) > 0) {
    split(line, a, /[|,]/)
    n++; px[n]=a[2]; py[n]=a[3]; pz[n]=a[4]
  }
}
{
  x0=$2; y0=$3; z0=$4
  mind=1e9
  for(i=1;i<=n;i++){
    dx=px[i]-x0; dy=py[i]-y0; dz=pz[i]-z0
    d=sqrt(dx*dx+dy*dy+dz*dz)
    if(d<mind) mind=d
  }
  printf "%s %.5f\n", $1, mind
}' "$targets"
