#!/usr/bin/perl

#  an ldac formatted file giving test data
$LDAC = shift();
#  an output of "hac -R", i.e., "DATASTEM.topk" 
#  must have exactly the same number of lines as $LDAC and docs corresponds
$RES = shift();
#  only look at this many results from $RES
$MAXK = shift();

if ( !defined($MAXK) ) {
    print STDERR "Need 3rd argument, max results per doc in $RES\n";
    exit(0);
}

open(L,"<$LDAC") or die "Cannot open LDAC file\n";
open(R,"<$RES") or die "Cannot open result topk file\n";

$totprec = 0;
$totrecall = 0;
$cnt = 0;
$zeros = 0;
while ( defined($ldac=<L>) ) {
    $res=<R>;
    if ( !defined($res) ) {
	print STDERR "ran out of results\n";
	exit(1);
    }
    #  cleanup input
    chomp($res);
    chomp($ldac);
    $res =~ s/^\s+//;
    $ldac =~ s/^\s+//;
    
    my @res = split(/\s+/,$res);
    #   remove any parts beyond MAXK
    splice @res,$MAXK,$#res-$MAXK+1;
    my %res = ();
    foreach my $r ( @res ) {
	$res{$r} = 1;
    }
    my @ldac = split(/\s+/,$ldac);
    my %ldac = ();
    shift(@ldac);
    foreach my $l ( @ldac ) {
	$l =~ s/:.*//;
	$ldac{$l} = 1;
    }
    # print "$cnt: l=$#ldac r=$#res\n";
    # print "RES: #" . join("#",@res) . "#\n";
    # print "LDAC: #" . join("#",@ldac) . "#\n";
    my $recall = 0;
    my $prec = 0;
    # print "PREC got: ";
    foreach my $r ( keys(%res) ) {
	if ( defined($ldac{$r}) ) {
	    # print " '$r'";
	    $prec ++;
	}
    }
    # print "\n";
    foreach my $l ( keys(%ldac) ) {
	if ( defined($res{$l}) ) {
	    $recall ++;
	}
    }
    # print "$cnt: r=$recall p=$prec\n";
    if ( $#ldac<0 ) {
	$zeros++;
	next;
    }
    my $restot = ($#ldac + 1.000000000001);
    if ( $restot > $MAXK) {
	$restot = $MAXK + 0.000000000001;
    }
    $recall /= ($#ldac + 1.00);
    $prec /= $restot;
    # print "$cnt: $recall (/$#ldac) $prec (/$restot)\n";
    $totrecall += $recall;
    $totprec += $prec;
    $cnt++;
}
close(L);
close(R);
$totrecall /= $cnt;
$totprec /= $cnt;

print "Zero test records = $zeros/$cnt\n";
print "Precision = $totprec, Recall = $totrecall\n";
