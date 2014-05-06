#!/usr/bin/perl

use utf8;
use POSIX;

$ALPHA = 4;
$POW = 1;

$cntfile = shift();
$stem = shift();
$startcnt = shift() - 1;

open(S,"<$cntfile");
while ( ($s=<S>) ) {
	my @a = split(/ /,$s);
	$wcs{lc($a[4])} = $a[3];
}
close(S);

open(R,"<$stem.tokens");
$line = 0;
$ptot = 0;
while ( defined($r=<R>) ) {
    chomp();
    $line++;
    if ( !defined($wcs{lc($r)}) ) {
	print STDERR "Cannot find word '$r'/$line in '$stem.tokens' \n";
	exit(1);
    }
    my $v = POSIX::pow($wcs{lc($r)} + $ALPHA, $POW);
    $pp[$line] = $v;
    $ptot += $v;
    # print "$v\n";
}
close(R);

for (my $i=1; $i<=$line; $i++) {
    $pp[$i] /= $ptot;
}

open(D,"<$stem.dit");
open(W,"<$stem.wit");
$ptot = 0;
$pcnt = 0;
while ( defined(my $d=<D>) ) {
    my $w=<W>;
    chomp($d);
    chomp($w);
    if ( $d>=$startcnt ) {
	# print "adding $w/$d of $pp[$w]\n";
	$ptot -= log($pp[$w]) / log(2.0);
	$pcnt++;
    }
}
close(W);
close(D);
$ptot /= $pcnt;
print "Perplexity = $ptot\n";
