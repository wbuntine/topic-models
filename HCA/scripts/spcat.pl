#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use POSIX;
use Getopt::Long;
use Pod::Usage;
           
GetOptions(
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 1 );

my $file1 = shift();
my $file2 = shift();

open(F1,"<$file1") or die "Cannot open '$file1': $!\n";
open(F2,"<$file2") or die "Cannot open '$file2': $!\n";

my $doc1 = int(<F1>);
my $doc2 = int(<F2>);
my $W1 = int(<F1>);
my $W2 = int(<F2>);
if ( $W1 != $W2 ) {
    print STDERR "Sparse files have different word counts\n";
    exit(1);
}
my $C1 = int(<F1>);
my $C2 = int(<F2>);

my $doc = $doc1+$doc2;
my $C = $C1+$C2;
print "$doc\n$W1\n$C\n";
while ( defined($_=<F1>) ) {
    print;
}
close(F1);
while ( defined($_=<F2>) ) {
    chomp();
    my @a = split();
    $a[0] += $doc1;
    print join(" ",@a) . "\n";
}
close(F2);

__END__

=head1 NAME
  
spcat.pl - merge two sparse matrix (docword) data files

=head1 SYNOPSIS
    
spcat.pl SFILE1 SFILE2 > MFILE

Options:

    SFILE1              first sparse matrix file
    SFILE2              second sparse matrix file
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Data format is sparse matrices ala
I<http://archive.ics.uci.edu/ml/datasets/Bag+of+Words>.
The two files must have the same vocabularies.  This merges the indices
for documents.  Output to standard out.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
