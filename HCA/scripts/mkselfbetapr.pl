#!/usr/bin/perl

#   use this to build a .betapr file from the data itself
#         sort -n   test.wit | uniq -c | ./mkselfbetapr.pl > test.betapr

use POSIX;

$ALPHA = 4.0;
$POW = 1.00;
$zero = POSIX::pow($ALPHA, $POW);

$last = 1;
while ( <> ) {
    @a = split();
    # print "$a[1] $a[0]\n";
    for (  ; $last<$a[1]; $last++) {
	print "$zero\n";
    }
    print POSIX::pow($a[0]+$ALPHA, $POW) . "\n";
    $last++;
}
