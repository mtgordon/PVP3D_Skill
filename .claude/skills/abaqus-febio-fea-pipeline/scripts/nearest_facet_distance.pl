#!/usr/bin/env perl
# For every point in targets.txt, find the minimum distance to the nearest
# POINT ON A FACET (triangle) in reference_facets.txt -- not just the
# nearest vertex. See references/geometry-cross-referencing.md Rule 8 for
# why this matters: a proximity-based facet contact (e.g. FEBio's
# tied-node-on-facet) projects onto a facet's surface, which can be
# meaningfully closer than any of that facet's own vertices, so a
# nearest-VERTEX measurement (nearest_point_distances.sh) can overstate the
# real contact gap by a wide margin.
#
# Usage: nearest_facet_distance.pl targets.txt nodes.txt facets.txt
#   targets.txt : "id x y z" or "id|x,y,z" -- points to measure from
#   nodes.txt   : "id x y z" or "id|x,y,z" -- node lookup table for the
#                 facets (build_node_lookup.sh output works directly)
#   facets.txt  : one facet per line, whitespace-separated node IDs
#                 (3 = triangle; 4 = quad, auto-split into 2 triangles)
#
# Output: "<target_id> <min_distance> <closest_facet_line_number>" per line

use strict; use warnings;

my ($targets_f, $nodes_f, $facets_f) = @ARGV;
die "usage: nearest_facet_distance.pl targets.txt nodes.txt facets.txt\n" unless $facets_f;

sub read_points {
    my ($f) = @_;
    my %pts;
    open(my $fh, '<', $f) or die "can't open $f: $!\n";
    while (<$fh>) {
        chomp;
        next unless /\S/;
        my ($id, $rest) = split(/[|\s]+/, $_, 2);
        my @xyz = split(/[,\s]+/, $rest);
        $pts{$id} = [@xyz[0..2]];
    }
    close $fh;
    return \%pts;
}

my $targets = read_points($targets_f);
my $nodes = read_points($nodes_f);

my @tris; # each: [n1,n2,n3] (node IDs)
open(my $ff, '<', $facets_f) or die "can't open $facets_f: $!\n";
while (<$ff>) {
    my @n = split;
    next unless @n >= 3;
    push @tris, [$n[0], $n[1], $n[2]];
    push @tris, [$n[0], $n[2], $n[3]] if @n >= 4;  # split quad into 2 tris
}
close $ff;

# closest point on triangle ABC to point P (Ericson, "Real-Time Collision
# Detection" 5.1.5) -- returns squared distance
sub sqdist_point_tri {
    my ($p, $a, $b, $c) = @_;
    my @ab = map { $b->[$_]-$a->[$_] } 0..2;
    my @ac = map { $c->[$_]-$a->[$_] } 0..2;
    my @ap = map { $p->[$_]-$a->[$_] } 0..2;
    my $d1 = $ab[0]*$ap[0]+$ab[1]*$ap[1]+$ab[2]*$ap[2];
    my $d2 = $ac[0]*$ap[0]+$ac[1]*$ap[1]+$ac[2]*$ap[2];
    if ($d1<=0 && $d2<=0) { return sqd($p,$a); }
    my @bp = map { $p->[$_]-$b->[$_] } 0..2;
    my $d3 = $ab[0]*$bp[0]+$ab[1]*$bp[1]+$ab[2]*$bp[2];
    my $d4 = $ac[0]*$bp[0]+$ac[1]*$bp[1]+$ac[2]*$bp[2];
    if ($d3>=0 && $d4<=$d3) { return sqd($p,$b); }
    my $vc = $d1*$d4-$d3*$d2;
    if ($vc<=0 && $d1>=0 && $d3<=0) {
        my $v = $d1/($d1-$d3);
        my @q = map { $a->[$_]+$v*$ab[$_] } 0..2;
        return sqd($p,\@q);
    }
    my @cp = map { $p->[$_]-$c->[$_] } 0..2;
    my $d5 = $ab[0]*$cp[0]+$ab[1]*$cp[1]+$ab[2]*$cp[2];
    my $d6 = $ac[0]*$cp[0]+$ac[1]*$cp[1]+$ac[2]*$cp[2];
    if ($d6>=0 && $d5<=$d6) { return sqd($p,$c); }
    my $vb = $d5*$d2-$d1*$d6;
    if ($vb<=0 && $d2>=0 && $d6<=0) {
        my $w = $d2/($d2-$d6);
        my @q = map { $a->[$_]+$w*$ac[$_] } 0..2;
        return sqd($p,\@q);
    }
    my $va = $d3*$d6-$d5*$d4;
    if ($va<=0 && ($d4-$d3)>=0 && ($d5-$d6)>=0) {
        my $w = ($d4-$d3)/(($d4-$d3)+($d5-$d6));
        my @q = map { $b->[$_]+$w*($c->[$_]-$b->[$_]) } 0..2;
        return sqd($p,\@q);
    }
    my $denom = 1.0/($va+$vb+$vc);
    my $v = $vb*$denom; my $w = $vc*$denom;
    my @q = map { $a->[$_]+$ab[$_]*$v+$ac[$_]*$w } 0..2;
    return sqd($p,\@q);
}
sub sqd { my ($p,$q)=@_; ($p->[0]-$q->[0])**2+($p->[1]-$q->[1])**2+($p->[2]-$q->[2])**2; }

for my $tid (sort keys %$targets) {
    my $p = $targets->{$tid};
    my ($best, $best_i);
    for my $i (0..$#tris) {
        my ($n1,$n2,$n3) = @{$tris[$i]};
        next unless exists $nodes->{$n1} && exists $nodes->{$n2} && exists $nodes->{$n3};
        my $d = sqdist_point_tri($p, $nodes->{$n1}, $nodes->{$n2}, $nodes->{$n3});
        if (!defined($best) || $d < $best) { $best = $d; $best_i = $i; }
    }
    printf "%s %.5f %d\n", $tid, sqrt($best // -1), $best_i // -1;
}
