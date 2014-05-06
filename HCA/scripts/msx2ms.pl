#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use utf8;
use POSIX;
use Getopt::Long;
use Pod::Usage;

# encoding pragmas follow any includes like "use"
use encoding 'utf8';
use open ':utf8';
binmode STDIN, ":utf8";
binmode STDERR, ":utf8";

#  by default use unnormalised PMI, ... set this for NPMI
my $usenorm = 0;

GetOptions(
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input cooc file")
      if ( $#ARGV != 0 );

my $stem = shift();

my %wind = ();
print STDERR "Reading word file\n";
open(C,"<${stem}.words") or die "Cannot open input '${stem}.words': $!";
while ( ($_=<C>) ) {
    if ( /^([0-9]+)\s[^\s]+\s[0-9a-fA-F]+\s[0-9]+\s[0-9]+\s(.*)$/ ) {
	chomp($2);
	$wind{$2} = $1;
    } else {
	die "Cannot parse line in '${stem}.words': $_";
    }
}
close(C);

open(S,">$stem.ms") or die "Cannot open output '$stem.ms': $!";
foreach my $TYPE ( ("n","v","j","r") ) {
    open(C,"<$stem.ms$TYPE") or die "Cannot open input '$stem.ms$TYPE': $!";
    while ( ($_=<C>) ) {
	if ( /^([^\s]+)\s([^\s]+)\s([^\s]+)$/ ) {
	    my $k1 = $wind{$1};
	    my $k2 = $wind{$2};
	    if ( defined($k1) && defined($k2) && $3 ne "0") {
		print S "$k1 $k2 $3\n";
	    }
	} else {
	    die "Cannot parse line in '$stem.words': $_";
	}
    }
    close(C);
}
close(S);

__END__

=head1 NAME
  
msx2ms.pl

=head1 SYNOPSIS
    
msx2ms.pl STEM

Options:

    STEM                Stem of set to get word counts.
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Read the similarity measures from the STEM.ms? files
("msr" for adverb, "msj" and adjective, "msn" for noun and "msv"
for verb).
This are listed by word.  Convert to
dictionary index form using STEM.words to get
indices, and represent as a sparse matrix.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
