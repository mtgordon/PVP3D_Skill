#!/usr/bin/env bash
# Check a FEBio .feb file for the node-ID-monotonicity gotcha (see
# references/febio-xml-format.md gotcha 1) BEFORE running febio4 on it.
#
# FEBio requires <node id="N"> declarations to appear in non-decreasing
# order of N across the entire <Mesh> section, regardless of which
# <Nodes name="..."> block they're in. Violating this fails with a
# confusing "invalid value for attribute id" error far from the real
# cause. This script flags every point where the ID sequence decreases,
# and which two <Nodes> blocks straddle the drop.
#
# Usage: check_monotonic_node_ids.sh model.feb

set -euo pipefail
f="$1"

awk '
  /<Nodes name=/ { blockname=$0; gsub(/.*name="/,"",blockname); gsub(/".*/,"",blockname) }
  match($0, /<node id="[0-9]+">/) {
    id = $0
    gsub(/.*<node id="/,"",id); gsub(/".*/,"",id); id+=0
    if (prev != "" && id < prev) {
      print "DROP at line " NR ": id " prev " (in \"" prevblock "\") -> id " id " (in \"" blockname "\")"
    }
    prev = id
    prevblock = blockname
  }
' "$f"

echo "Done. No output above lines means node IDs are monotonically non-decreasing."
