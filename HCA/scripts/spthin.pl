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
      if ( $#ARGV != 3 );

my $file = shift();
my $tokens = shift();
my $words = shift();
my $stem = shift();

open(F,"<$file") or die "Cannot open '$file': $!\n";

my $doc = int(<F>);
my $W = int(<F>);
my $C = int(<F>);

if ( $words>=$W+10 ) {
	print STDERR "Not enough words to cut out $words/$W\n";
	exit(1);
}

#   lines this word has
my @lines = ();
#   count this word has
my @count = ();
my @vocab = ();

for (my $i=0; $i<$W; $i++) { $count[$i]=0; $lines[$i]=0; }

#   first pass to count the vocab
while ( defined($_=<F>) ) {
	chomp();
	my @a = split();
	$lines[$a[1]]++;
	$count[$a[1]] += $a[2];
}
close(F);

open(F,"<$tokens") or die "Cannot open '$tokens': $!\n";
my $ll = 0;
while ( defined($_=<F>) ) {
    chomp();
    $vocab[$ll] = $_;
    $ll++;
}
close(F);
    
my @wordsort = sort { $count[$b] <=> $count[$a] } (1..($W));

my $totlines = 0;
my %wordset = ();
my @wordmap = ();

open(F,">$stem.vocab") or die "Cannot open '$stem.vocab': $!\n";
for (my $i=0; $i<$words; $i++) {
    $totlines += $lines[$wordsort[$i]];
    print F $vocab[$wordsort[$i]] . "\n";
    $wordmap[$wordsort[$i]] = $i+1;
}
close(F);

open(T,">$stem.docword") or die "Cannot open '$stem.docword': $!\n";
open(F,"<$file") or die "Cannot open '$file': $!\n";
<F>; <F>; <F>;

print T "$doc\n$words\n$totlines\n";
my $lines = 0;
while ( defined($_=<F>) ) {
    chomp();
    my @a = split();
    if ( defined($wordmap[$a[1]]) ) {
	print T "$a[0] $wordmap[$a[1]] $a[2]\n";
	$lines++;
    }
}
close(F);
close(T);
if ( $lines != $totlines ) {
    print STDERR "Line counts don't agree\n";
}

__END__

=head1 NAME
  
spthin.pl - thin out vocabulary of sparse matrix (docword) file

=head1 SYNOPSIS
    
spthin.pl SFILE TFILE N STEM

Options:

    SFILE               sparse matrix file
    TFILE               vocab file
    N			number of most frequent words to take
    STEM                for output of STEM.docword, STEM.vocab
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Data format is sparse matrices ala
I<http://archive.ics.uci.edu/ml/datasets/Bag+of+Words>.
Thin out vocabulary from a sparse matrix file by taking top N words.
Output to other files.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
