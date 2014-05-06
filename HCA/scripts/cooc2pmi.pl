#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use utf8;
use POSIX;
use Getopt::Long;
use Pod::Usage;

my $ALPHA = 1.0;
my $donorm = 0;
my $dowin = 0;
my $SIGMA = 0;

# encoding pragmas follow any includes like "use"
use encoding 'utf8';
use open ':utf8';
binmode STDIN, ":utf8";
binmode STDERR, ":utf8";

#  by default use unnormalised PMI, ... set this for NPMI
my $usenorm = 0;

GetOptions(
    'alpha=s'     => \$ALPHA,
    'sigma=s'     => sub {
	local *_ = \$_[1];
	  /^([0-9\.]+)$/
	       or die("Invalid format for option -sigma.\n");
	   $SIGMA = $1;
    },
    'norm!'      =>  \$donorm,
    'win!'       =>  \$dowin,
    'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
    'norm!'        => \$usenorm,
    'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input cooc file")
      if ( $#ARGV != 1 );

if ( $SIGMA>0 ) {
    $ALPHA = 0;
    $donorm = 0;
}

#   the total is in '-1'
my $coocfile = shift();
my $stem = shift();

my $fullcnt = 0;
my $fullwin = 0;
my @totcnt = ();   # from cooc file
my @totwin = ();   # from cooc file
my @cnt = ();      # from .words file
my $cnttot = 0;

#  first pass fills totals by matching k2==-1
print STDERR "Reading $coocfile for totals\n";
open(C,"<$coocfile") or die "Cannot open input '$coocfile': $!";
$_ = <C>;
if ( /^-1\s-1\s([0-9\.]+)\s([0-9\.]+)/ ) {
    $fullcnt = $1;
    $fullwin = $2;
} else {
    print STDERR "Cannot parse totals line for '$coocfile'\n";
    exit(1);
}
while ( ($_=<C>) ) {
    if ( /^([0-9]+)\s+-1\s([0-9\.]+)\s([0-9\.]+)/ ) {
	if ( $dowin ) {
	    $cnt[$1] = $3+$ALPHA;
	} else {
	    $cnt[$1] = $2+$ALPHA;
	}
	$cnttot += $cnt[$1];
    }
}
close(C);

print STDERR "Using cnttot=$cnttot versus fullcnt=$fullcnt\n";

print STDERR "Reading $coocfile for PMI\n";

open(C,"<$coocfile") or die "Cannot reopen input '$coocfile': $!";
while ( ($_=<C>) ) {
    if ( /^([0-9]+)\s+([0-9]+)\s+([0-9\.]+)/ ) {
	my $i1 = $1;
	my $i2 = $2;
	my $if = $3;
	if ( defined($cnt[$i1]) && defined($cnt[$i2]) ) {
	    my $val = $if / $fullwin;
	    my $norm = $val;
	    $val /= ($cnt[$i1]/$cnttot)*($cnt[$i2]/$cnttot);
	    if ( $SIGMA>0) {
		#  shrink using +/- $SIGMA*sqrt() back to 1
		if ( $val > 1 ) {
		    $val = ($if - $SIGMA*sqrt($if)) / $fullwin;
		    $val /= (($cnt[$i1]+$SIGMA*sqrt($cnt[$i1]))/$cnttot)
			* (($cnt[$i2]+$SIGMA*sqrt($cnt[$i2]))/$cnttot);
		    if ( $val< 1 ) {
			$val = 1;
		    } 
		} else {
		    $val = ($if + $SIGMA*sqrt($if)) / $fullwin;
		    $val /= (($cnt[$i1]-$SIGMA*sqrt($cnt[$i1]))/$cnttot)
			* (($cnt[$i2]-$SIGMA*sqrt($cnt[$i2]))/$cnttot);
		    if ( $val> 1 ) {
			$val = 1;
		    } 
		}
	    }
	    $val = log($val);
	    if ( $donorm ) {
		$val /= -log($norm);
	    }
	    if ( $val != 1 ) {
		printf "$i1 $i2 %.4f\n", $val;
	    }
	}
    }
}
close(C);

__END__

=head1 NAME
  
cooc2pmi.pl - convert ".cooc" file to ".pmi" to stdout.

=head1 SYNOPSIS
    
cooc2pmi.pl [--alpha A --sigma S --norm --win] COOCFILE STEM > PMIFILE

Options:

    COOCFILE            Source co-occurrence counts file, -1 for total cols
    STEM                Stem of set to get word counts.
    --alpha A           Alpha offset for counts, default is 4.0
    -h, --help          display help message and exit.
    --man              print man page and exit.
    --norm             use the newer normalised PMI score
    --sigma S          shrink counts back to independence by this many
    --win              for single word probs, use window probability,
                        not the default occurence probability

=head1 DESCRIPTION

Convert co-occurrence count file into a PMI file.   Must have the total
entries (lines "W -1 Cnt").  Note count matrices have word indices
offset by 0.  The co-occurrence count file should have been generated
by
I<linkCoco> from the DCA-Bags release.

=head1 SEE ALSO

I<linkCoco>(1).

DCA-bags and ddca websites are inside
F<http://forge.nicta.com.au>

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011-2013 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
