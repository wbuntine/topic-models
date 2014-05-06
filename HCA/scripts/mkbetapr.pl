#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use utf8;
use POSIX;
use Getopt::Long;
use Pod::Usage;

my $ALPHA = 4.0;
my $POW = 1.0;
my $zero = POSIX::pow($ALPHA, $POW);
# set this to subtract own counts from data
my $subtract = 0;

# encoding pragmas follow any includes like "use"
use encoding 'utf8';
use open ':utf8';
binmode STDIN, ":utf8";
binmode STDERR, ":utf8";

GetOptions(
      'alpha=s'     => \$ALPHA,
      'pow=s'     => \$POW,
      'subtract!'   => \$subtract,
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 1 );

# .words file to use as counts data
my $mstfile = shift();
my $stem = shift();

#  hash table of counts for tokens
my %wcs = ();

print STDERR "Using '$mstfile' for statistics and matching tokens\n";
open(S,"<$mstfile");
while ( (my $s=<S>) ) {
    chomp($s);
    if ( $s =~ /^[0-9]+ text [^ ]+ ([0-9]+) (.*)/ ) {
	$wcs{lc($2)} = $1;
    }
}
close(S);

if ( -f "$stem.words" ) {
    print STDERR "Using '$stem.words' for own tokens and statistics\n";
    open(R,"<$stem.words");
    while ( (my $s=<R>) ) {
	chomp($s);
	if ( $s =~ /^[0-9]+ text [^ ]+ ([0-9]+) (.*)/ ) {
	    my $cntm = $wcs{lc($2)};
	    my $cnt = $1;
	    if ( !defined($cntm) ) {
		print STDERR "Cannot find word '$2' in file '$mstfile' \n";
		$cntm = 0;
	    }
	    if ( $subtract ) {
		$cntm -= $cnt;
		if ( $cntm<0 ) { $cntm = 0; }
	    }
	    my $v = POSIX::pow($cntm + $ALPHA, $POW);
	    print "$v\n";
	}
    }
} else {
    if ( $subtract ) {
	print STDERR "No '$stem.words' file so cannot subtract\n";
	exit(1);
    }
    print STDERR "Using '$stem.tokens' for own tokens\n";
    open(R,"<$stem.tokens") or die "Cannot open tokens file\n";
    while ( ($_=<R>) ) {
	chomp();
	my $cntm = $wcs{lc($_)};
	if ( !defined($cntm) ) {
	    print STDERR "Cannot find word '$_' in file '$mstfile' \n";
	    $cntm = 0;
	}
	my $v = POSIX::pow($cntm + $ALPHA, $POW);
	print "$v\n";
    }
}
close(R);


__END__

=head1 NAME
  
mkbetapr.pl - convert ".word" files to a beta prior to stdout

=head1 SYNOPSIS
    
mkbetapr.pl [--alpha A|--pow P] SRCWORD STEM

Options:

    SRCWORD             Source ".words" file giving counts and tokens
    STEM                Stem for self
    --alpha A           Alpha offset for counts, default is 4.0
    --pow P             Raise to power to dampen, default is 1.0
    --subtract          Subtract own counts if source was superset of data
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Creates a informed weights for the beta prior.  Must match tokens
so use a source "STEM.tokens" to arrange matching.  Uses a source collections
".words" file to get details.  With the --subtract option must have "STEM.words"
file.  Output counts to standard out.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
