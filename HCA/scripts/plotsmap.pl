#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell
#
#   usage:
#         plotsmap.pl  7,23,123 60 < t1.smap 

use POSIX;
use Getopt::Long;
use Pod::Usage;
           
GetOptions(
      'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 1 );

@w = split(/,/,shift());
$K = shift();
foreach my $w ( @w ) {
    $w{$w} = 1;
}
$n_w = $#w+1;

$maxv = 0;
$w = 0;
while ( <> ) {
    if ( /^([^\(]+)\(([0-9]+)\)\: ([0-9\.\/ ]+)perp=([0-9\.]+)/ 
	 && defined($w{$2}) ) {
	$name[$w] = $1;
	$ind[$w] = $2;
	$perp[$w] = $4;
	@a = split(/ /,$3);
	$norm = 0;
	for (my $i=0; $i<=$#a; $i++) {
	    ($k,$v) = split(/\//,$a[$i]);
	    $norm += $v;
	}
	$n = 0;
	open(F,">$ind[$w].gpd");
	for (my $i=0; $i<=$#a; $i++) {
	    ($k,$v) = split(/\//,$a[$i]);
	    for ( ; $n<$k; $n++) {
		print F "$n 0\n";
	    }
	    $n = $k+1;
	    $v /= $norm;
	    if ( $v >= $maxv ) { $maxv = $v; }
	    print F "$k $v\n";
	}
	for ( ; $n<$K; $n++) {
	    print F "$n 0\n";
	}
	close(F);
	print STDERR "Word '$name[$w]' ($ind[$w]) ent=$perp[$w]\n";
	$w++;
    }
}
$scale = sprintf("%.1f", $maxv);

open(T,">try.gp");
print T "set terminal jpeg large\n";
print T "set output \"try.jpeg\"\n";
print T "set boxwidth 0.9 relative\n";
print T "set style fill solid 1.0\n";
print T "set multiplot layout 1,$w\n";
print T "set lmargin 0\n";
print T "set rmargin 0\n";
print T "set tmargin 3\n";
print T "set bmargin 3\n";
print T "set yrange [0:$scale]\n";
print T "unset xtics\n";
print T "unset ytics\n";
#
for (my $i=0; $i<$w; $i++) {
    print T "set title \"$name[$i]($ind[$i]) p=$perp[$i]\"\n";
    print T "unset key\n";
    print T "plot \"$ind[$i].gpd\" with boxes\n";
    
}
print T "unset multiplot\n";
close(T);
system("gnuplot < try.gp")==0 or die "Cannot run gnuplot\n";
print "Plot done as 'try.jpeg'\n";

for (my $i=0; $i<$w; $i++) {
	unlink("$ind[$i].gpd");
}
unlink("try.gp");

__END__

=head1 NAME
  
plotsmap.pl - use gnuplot to plot topic probabilities for a list of words

=head1 SYNOPSIS
    
plotsmap.pl WORDLIST TOPICCOUNT < RESSTEM.smap


Options:

    WORDLIST            Comma separated list of word indices, no spaces
    TOPICCOUNT          K, or total number of topics
    RESSTEM.smap        output ".smap" file from running hca with '-lsp,N,M'
    -h, --help          display help message and exit.
     --man              print man page and exit.

=head1 DESCRIPTION

Creates a jpeg file
F<try.jpeg> that plots the word probabilities
for the words in the list.  Since they are on one plot, best not use
more than 4 words at once.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
