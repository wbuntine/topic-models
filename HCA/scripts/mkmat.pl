#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use POSIX;
use Getopt::Long;
use Pod::Usage;
           
# encoding pragmas follow any includes like "use"
use encoding 'utf8';
use open ':utf8';
binmode STDIN, ":utf8";
binmode STDERR, ":utf8";

my $prec = 3;

GetOptions(
      'prec=i'   => \$prec,
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 2 );

my $vocfile = shift();  # vocab. for source matrix
my $mtxfile = shift();  # source matrix (PMI of COOC)
my $stem = shift();     #  output stem, but also stem for desired vocab

my %voc = ();
my @map = ();

open(T,"<$stem.tokens");
my $cnt = 0;
while ( ($_=<T>) ) {
    chomp();
    $voc{$_} = $cnt;
    $cnt++;
}
close(T);

print STDERR "Read $cnt vocab from '$stem.tokens'\n";

#  build map from source vocab

open(T,"<$vocfile");
$cnt = 0;
my $mapped = 0;
my %gotv = ();
while ( ($_=<T>) ) {
    chomp();
    if ( defined($voc{$_}) ) {
	if ( !defined($gotv{$voc{$_}}) ) {
	    $map[$cnt] = $voc{$_};
	    $gotv{$voc{$_}} = 1;
            $mapped++;
	}
    } else {
	$map[$cnt] = -1;
    }
    $cnt++;
}
close(T);

print STDERR "Read $cnt vocab from '$vocfile'\n";
print STDERR "Mapped $mapped vocab from '$vocfile' to '$stem'\n";

open(M,"<$mtxfile");
while ( ($_=<M>) ) {
    chomp();
    if ( /^([^ ]+) ([^ ]+) ([^ ]+)$/ ) {
	if ( !defined($map[$1]) ) {
	    print STDERR "Map undefined for 1:$1\n";
        } elsif ( !defined($map[$2]) ) {
	    print STDERR "Map undefined for 2:$2\n";
        } elsif ( $map[$1]>=0 && $map[$2]>=0 ) {
            my $score = sprintf "%.${prec}f", $3;
	    print "$map[$1] $map[$2] $score\n";
	}
    } else {
	print STDERR "Bad line: $_\n";
	exit(1);
    }
}
close(M);

__END__

=head1 NAME
  
mkmat.pl - convert a similarity matrix from one 

=head1 SYNOPSIS
    
mkmat.pl SRCVOC SRCMTX STEM

Options:

    SRCVOC              Source vocab, one token per line
    SRCMTX              Matrix in sparse format to convert
    STEM                Stem for output file and matching token file
    -h, --help          display help message and exit.
     --man              print man page and exit.
    --prec P            decimal places after '.' to keep

=head1 DESCRIPTION

Convert the matrix in the file SRCMTX to be suitable for the
data set with stem STEM.  SRCMTX has matching vocabulary in the token file
SRCVOC.  Reads STEM.tokens to convert the matrix.   Output to stdout.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
