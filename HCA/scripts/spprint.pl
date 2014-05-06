#!/usr/bin/perl -w

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

GetOptions(
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 1 );

my $vocfile = shift();
my $mtxfile = shift();

my @voc = ();

open(T,"<$vocfile");
my $cnt = 1;
while ( ($_=<T>) ) {
    chomp();
    $voc[$cnt] = $_;
    $cnt++;
}
close(T);

open(M,"<$mtxfile");
while ( ($_=<M>) ) {
    chomp();
    if ( /^([^ ]+) ([^ ]+) ([^ ]+)$/ ) {
	print "$voc[$1] $voc[$2] $3\n";
    } 
}
close(M);

__END__

=head1 NAME
  
prmat.pl - print a similarity matrix (given in sparse, docword, format)

=head1 SYNOPSIS
    
mkmat.pl VOC MTX

Options:

    VOC              Vocab, one token per line
    MTX              Matrix in sparse format to print
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Data format is sparse matrices ala
I<http://archive.ics.uci.edu/ml/datasets/Bag+of+Words>.
Print with tokens instead of indices.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
