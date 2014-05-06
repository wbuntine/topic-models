#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use List::Util;

my $ldac = 0;
my $seed = 0;
           
GetOptions(
      'seed=i'  => \$seed,
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 2 );

my $stem1 = shift();
my $stem2 = shift();
my $TRAIN = shift();
my $N = 0;

if ( $seed>0 ) {
    srand($seed);
}

if ( -e "$stem1.ldac" ) {
    $ldac = 1;
    if ( ! -e "$stem1.srcpar" ) {
	print STDERR "Cannot read '$stem1.srcpar'\n";
	exit(1);
    }
    open(L,"grep documents $stem1.srcpar |");
    $N = <L>;
    chomp($N);
    $N =~ s/.*=//;
    $N = int($N);
    close(L);
} elsif ( -e "$stem1.txtbag" ) {
    open(L,"head -1 $stem1.txtbag |");
    $N = int(<L>);
    close(L);
} else {
    print STDERR "Cannot locate bag file '$stem1.ldac' or '$stem1.txtbag'\n";
    exit(1);
}
if ( !defined($N) || $N<=0 ) {
    print STDERR "Cannot get document count\n";
    exit(1);
}
print STDERR "Extracting $TRAIN training docs from $N\n";

# read epoch file and store in array to by indexed by doc index
open(L,"<$stem1.epoch") or
    die "Cannot read '$stem1.epoch'\n";
my $E = int(<L>);
my @epoch = ();
for (my $e = 0; $e<$E; $e++) {
    my $n = int(<L>);
    for ( ; $n>0; $n--) {
	push(@epoch,$e);
    }
}
close(L);
if ( $#epoch+1 != $N ) {
    print STDERR "Doc count in epoch file doesnt match\n";
    exit(1);
}

# make a hash giving the training set
my @rlist = List::Util::shuffle(0 .. ($N-1));
my %train = ();
for (my $i=0; $i<$TRAIN; $i++) {
    $train{$rlist[$i]} = 1;
}
@rlist = ();

#  build training part

if ( $ldac ) {
    open(I,"<$stem1.ldac");
    open(O,">$stem2.ldac");
} else {
    open(I,"<$stem1.txtbag");
    open(O,">$stem2.txtbag");
    print O <I>;
    print O <I>;
}
open(E,">$stem2.epoch");

my @eout = ();
for (my $e = 0; $e<$E; $e++) {
    $eout[$e] = 0;
}
for (my $i=0; $i<$N; $i++) {
    my $line = <I>;
    if ( $train{$i} ) {
	$eout[$epoch[$i]]++;
	print O $line;
    }
}
print E "$E\n";
for (my $e = 0; $e<$E; $e++) {
    print E " $eout[$e]\n";
}
seek I,0,0;
@eout = ();
for (my $e = 0; $e<$E; $e++) {
    $eout[$e] = 0;
}
for (my $i=0; $i<$N; $i++) {
    my $line = <I>;
    if ( !defined($train{$i}) ) {
	$eout[$epoch[$i]]++;
	print O $line;
    }
}
print E "$E\n";
for (my $e = 0; $e<$E; $e++) {
    print E " $eout[$e]\n";
}
close(E);
close(I);
close(O);

__END__

=head1 NAME
  
berand.pl - build train/test of the bag file and create matching epoch file

=head1 SYNOPSIS
    
spcat.pl INSTEM OUTSTEM TRAIN

Options:

    INSTEM              file stem for input data set
    OUTSTEM             file stem for output data set, only subset of files
    TRAIN               size of training
    -h, --help          display help message and exit.
     --man              print man page and exit.
    --seed S            seed ranom number generator

=head1 DESCRIPTION

First the program looks for 
F<INSTEM.ldac> and then
F<INSTEM.txtbag>.
Which ever is found, then a training set of size TRAIN and testing set,
the remainder, are extracted and placed in
F<OUTSTEM.ldac> or
F<OUTSTEM.txtbag>.
Apart from the train test split, the data order is retained.
The epoch file is then converted to represent the split
and output to 
F<OUTSTEM.epoch>.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
