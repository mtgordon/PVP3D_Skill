#!/usr/bin/env bash
# Extract a flat "id|x,y,z" lookup table from a FEBio .feb file's
# <node id="N">x,y,z</node> lines, for fast coordinate cross-referencing
# (see references/geometry-cross-referencing.md Rules 2-3). Build this
# once per file, then do all subsequent nearest-point/exact-match queries
# against the resulting table with awk associative arrays instead of
# re-scanning the source file per query.
#
# Usage: build_node_lookup.sh model.feb > model_nodes.txt
#
# Output format: one line per node, "id|x,y,z"

set -euo pipefail
f="$1"

grep -oE '<node id="[0-9]+">[^<]+</node>' "$f" |
  sed -E 's/<node id="([0-9]+)">([^<]+)<\/node>/\1|\2/'
