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

my $file = shift();
my $docs = shift();

open(F,"<$file") or die "Cannot open '$file': $!\n";

my $doc = int(<F>);
my $W = int(<F>);
my $C = int(<F>);

if ( $doc<$docs ) {
	print STDERR "Not enough documents to cut out $docs\n";
	exit(1);
}

my $rd = 0;
my $lines = 0;
while ( defined($_=<F>) && $rd<=$docs ) {
	chomp();
	my @a = split();
	$rd = $a[0];
	$lines ++;
}
if ( $rd != $docs+1 ) {
    	print STDERR "Couldn't find place to cut out $docs\n";
	exit(1);
}
$lines -= 1;
seek(F,0,0);
<F>;
<F>;
<F>;
print "$docs\n$W\n$lines\n";
while ( defined($_=<F>) && $lines>0 ) {
    print;
    $lines--;
}
close(F);

__END__

=head1 NAME
  
spcut.pl - take first block of records out of sparse matrix (docword) file

=head1 SYNOPSIS
    
spcut.pl SFILE N 

Options:

    SFILE               sparse matrix file
    N			number of initial records to take
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Data format is sparse matrices ala
I<http://archive.ics.uci.edu/ml/datasets/Bag+of+Words>.
Take first N records out of sparse file.  Output to standard out.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
