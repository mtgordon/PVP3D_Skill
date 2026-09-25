#!/usr/bin/env perl
# Fit an N-term incompressible Ogden model to an Abaqus *Uniaxial Test Data
# table, via multi-start coordinate-descent least squares (no external
# numerical library needed -- useful when Python/scipy isn't available, or
# isn't working: on at least one Windows machine the default `python`/`py`
# resolved to an Anaconda-linked install with a BROKEN numpy (DLL load
# failure on `import numpy`), while a second installed version (`py -3.10`)
# had a perfectly good numpy/scipy -- if you ever do reach for Python here,
# check `py -0` for alternate installed versions before concluding
# Python/scipy just isn't usable at all.
#
# sigma_nominal(lambda) = Sum_i (2*mu_i/alpha_i) * (lambda^(alpha_i-1) - lambda^(-alpha_i/2-1))
# lambda = 1 + nominal_strain
#
# Usage: fit_ogden.pl data.txt [max_N]
#   data.txt : one "stress strain" pair per line (whitespace-separated),
#              matching Abaqus's *Uniaxial Test Data column order (stress
#              first). DOUBLE CHECK this against the actual source: a
#              plainly-listed "strain, stress" table (the more common way
#              data shows up outside an Abaqus *Uniaxial Test Data block)
#              has the OPPOSITE column order and will silently fit backwards
#              if fed in as-is -- the pair values alone don't disambiguate
#              which column is which, only knowing the source's convention
#              does. Include the (0,0) origin point if present; it doesn't
#              affect the fit either way.
#   max_N (optional, default 3): fits N=1..max_N and reports every one --
#              compare R^2, the stability flag, and the monotonicity flag
#              across them before picking; don't just take the largest
#              max_N by default (see the printed guidance at the end).
#
# Output, per N from 1 to max_N: fitted (mu_i, alpha_i) pairs (Abaqus/
# Ogden-1972 convention), the converted FEBio (c_i, m_i) (c_i=2*mu_i,
# m_i=alpha_i -- see febio-xml-format.md gotcha 16, including its bulk-
# modulus/Poisson's-ratio conversion for FEBio's required k), R^2, whether
# every term is Drucker-stable (mu_i*alpha_i > 0 for all i), whether the
# fitted curve is monotonic over the data range plus 10% beyond (a
# non-monotonic dip/overshoot is a stronger overfitting tell than R^2 alone
# -- R^2 can look great right up to the point a fit starts threading the
# data with a curve that dips or reverses direction between points), and a
# fit-vs-data comparison table.
#
# A concave-down/softening stress-strain shape (large initial slope,
# flattening out -- as opposed to a stiffening/hardening shape) commonly has
# NO fully term-wise-stable fit at any N: the unconstrained least-squares
# optimum needs alpha < 0 just to reach a decent R^2 even at N=1. That's not
# a bug in this script -- don't force artificially-restricted signs/ranges
# to chase stability on a shape that doesn't support it; report the
# instability honestly as a caveat on the extrapolation range instead (see
# the printed guidance, and febio-xml-format.md gotcha 16).

use strict; use warnings;

my $datafile = shift @ARGV or die "usage: fit_ogden.pl data.txt [max_N]\n";
my $max_n = shift @ARGV || 3;

open(my $fh, '<', $datafile) or die "can't open $datafile: $!\n";
my @data;
while (<$fh>) {
    next unless /(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)/;
    push @data, [$1+0, $2+0];  # [stress, strain]
}
close $fh;
die "need at least 2 data points\n" if @data < 2;

sub ogden_stress {
    # $params: flat arrayref mu_1,alpha_1,mu_2,alpha_2,...
    my ($params, $lambda) = @_;
    my $s = 0;
    for (my $i = 0; $i < @$params; $i += 2) {
        my ($mu, $alpha) = ($params->[$i], $params->[$i+1]);
        $s += (2*$mu/$alpha) * ($lambda**($alpha-1) - $lambda**(-$alpha/2-1));
    }
    return $s;
}

sub sse {
    my ($params) = @_;
    my $err = 0;
    for my $d (@data) {
        my ($stress, $strain) = @$d;
        my $pred = ogden_stress($params, 1 + $strain);
        $err += ($pred - $stress)**2;
    }
    return $err;
}

sub mean_stress { my $s = 0; $s += $_->[0] for @data; return $s / @data; }

sub r_squared {
    my ($params) = @_;
    my $ss_res = sse($params);
    my $mbar = mean_stress();
    my $ss_tot = 0;
    $ss_tot += ($_->[0] - $mbar)**2 for @data;
    return $ss_tot > 0 ? 1 - $ss_res / $ss_tot : 1;
}

# coordinate descent with shrinking step, from a given starting point
sub refine {
    my @params = @_;
    my $n = @params;
    my @step = map { $_ % 2 == 0 ? 0.02 : 1.0 } 0 .. $n - 1;  # mu step, alpha step
    for (1 .. 300) {
        my $improved = 0;
        my $cur = sse(\@params);
        for my $i (0 .. $n - 1) {
            for my $d (-1, 1) {
                my @trial = @params;
                $trial[$i] += $d * $step[$i];
                next if $i % 2 == 0 && abs($trial[$i]) < 1e-6;  # mu shouldn't cross 0
                next if $i % 2 == 1 && abs($trial[$i]) < 1e-3;  # nor alpha
                my $e = sse(\@trial);
                if ($e < $cur) { @params = @trial; $cur = $e; $improved = 1; }
            }
        }
        unless ($improved) { $_ /= 2 for @step; }
        last if $step[0] < 1e-9;
    }
    return @params;
}

srand(12345);  # deterministic across runs (this seeds rand(), unrelated to
               # Perl's separate hash-iteration randomization -- see
               # lessons-learned.md -- but pinned anyway so re-running the
               # same data gives the same fit)

for my $n (1 .. $max_n) {
    my ($best_params, $best_sse);
    my $n_restarts = 40 * $n;  # more terms -> more local optima -> more restarts
    for (1 .. $n_restarts) {
        my @p0;
        for (1 .. $n) {
            my $alpha0 = (rand() < 0.5 ? -1 : 1) * (0.3 + rand() * 15);
            my $mu0    = (rand() < 0.5 ? -1 : 1) * (0.05 + rand() * 5);
            push @p0, $mu0, $alpha0;
        }
        my @fit = refine(@p0);
        # coordinate descent alone tends to stall on a shallow ridge that
        # needs several parameters to move together to escape (empirically
        # ~2-3% worse R^2 than a gradient-based fit on a real dataset with
        # this script's very first version) -- a few rounds of "jiggle every
        # parameter at once, then re-refine, keep if better" recovers most
        # of that gap at a modest extra cost per restart
        for (1 .. 5) {
            my @jiggled = map { $fit[$_] + (rand() - 0.5) * 2 * ($_ % 2 == 0 ? 0.3 : 1.5) } 0 .. $#fit;
            my @repolished = refine(@jiggled);
            @fit = @repolished if sse(\@repolished) < sse(\@fit);
        }
        my $e = sse(\@fit);
        if (!defined($best_sse) || $e < $best_sse) { $best_sse = $e; $best_params = [@fit]; }
    }

    my $r2 = r_squared($best_params);
    my $stable = 1;
    for (my $i = 0; $i < @$best_params; $i += 2) {
        $stable = 0 if $best_params->[$i] * $best_params->[$i + 1] <= 0;
    }

    # monotonicity check over the data range plus 10% beyond
    my ($min_strain, $max_strain) = (1e9, -1e9);
    for my $d (@data) {
        $min_strain = $d->[1] if $d->[1] < $min_strain;
        $max_strain = $d->[1] if $d->[1] > $max_strain;
    }
    my $span = $max_strain - $min_strain;
    my $monotonic = 1;
    my $prev;
    for my $k (0 .. 60) {
        my $e = $min_strain + ($span * 1.1) * $k / 60;
        my $s = ogden_stress($best_params, 1 + $e);
        $monotonic = 0 if defined($prev) && $s < $prev - 1e-9;
        $prev = $s;
    }

    print "=" x 70, "\n";
    printf "N=%d term%s   R^2=%.6f   stable=%s   monotonic(data range +10%%)=%s\n",
        $n, $n == 1 ? "" : "s", $r2, $stable ? "yes" : "NO",
        $monotonic ? "yes" : "NO (overfitting signal -- see closing note)";
    for (my $i = 0; $i < @$best_params; $i += 2) {
        my ($mu, $alpha) = ($best_params->[$i], $best_params->[$i + 1]);
        printf "  term %d:  mu=%.6f  alpha=%.4f  (mu*alpha=%.3f)  ->  FEBio c%d=%.6f m%d=%.4f\n",
            $i / 2 + 1, $mu, $alpha, $mu * $alpha, $i / 2 + 1, 2 * $mu, $i / 2 + 1, $alpha;
    }
    print "  fit vs data:\n";
    for my $d (@data) {
        my ($stress, $strain) = @$d;
        my $pred = ogden_stress($best_params, 1 + $strain);
        printf "    strain=%.4f  data=%.5f  fit=%.5f  (%.1f%% off)\n",
            $strain, $stress, $pred, $stress != 0 ? 100 * abs($pred - $stress) / abs($stress) : 0;
    }
}

print "=" x 70, "\n";
print <<'NOTE';
Pick the smallest N where R^2 has saturated (further N buys little R^2 but a
jump in parameter magnitude/instability) over the N that fits best in
isolation -- with only a handful of data points, more terms than the data
supports overfits easily. The monotonic flag is a stronger overfitting tell
than R^2 alone: R^2 can look excellent right up to a fit that dips or
reverses direction between data points, which the R^2 number itself doesn't
show but a plotted curve does.

Term-wise instability (mu_i*alpha_i <= 0 for some i) is common and not
automatically disqualifying on its own -- but a genuinely concave-down/
softening stress-strain shape may have NO fully stable fit at any N. Don't
force a stable-looking fit that badly underfits just to get stability;
document the instability as a known extrapolation caveat instead (don't
trust the fit far outside the strain range the source data actually spans --
see febio-xml-format.md gotcha 16).
NOTE
