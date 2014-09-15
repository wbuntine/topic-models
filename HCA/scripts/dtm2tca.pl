#!/usr/bin/perl

$DTMDIR = shift();
$TCASTEM = shift();

if ( ! -e "$DTMDIR/lda-seq/info.dat" ) {
    print STDERR "Cmd arg #1 should be directory with entry 'lda-seq/info.dat'\n";
    exit(1);
}
if ( ! -e "$TCASTEM.par" ) {
    print STDERR "Cmd arg #2 should be a valid already-run 'tca' stem\n";
    exit(1);
}

open(I,"<$DTMDIR/lda-seq/info.dat") or die "Cannot open info.dat\n";

my $conc = 0;
my @alpha = ();
my $E = 0;
my $K = 0;
my $W = 0;
while ( ($_=<I>) ) {
    if ( /^ALPHA\s/ ) {
	chomp();
	my @a = split(/\s+/,$_);
	shift(@a);
	$K = shift(@a);
	if ( $K != $#a+1 ) {
	    print STDERR "Dimensions for ALPHA wrong in 'info.dat'\n";
	    exit(1);
	}
	my $tot = 0;
	for (my $i=0; $i<=$#a; $i++) {
	    $tot += $a[$i];
	}
	for (my $i=0; $i<=$#a; $i++) {
	    $a[$i] /= $tot;
	}
	#  $tot = conc. for the DP, @a is the base distribution
	$conc = $tot;
	@alpha = @a;
	print STDERR "Read alpha: K=$K, conc=$tot\n";
    }
    if ( /^SEQ_LENGTH\s/ ) {
	chomp();
	$E = $_;
	$E =~ s/^SEQ_LENGTH\s+//;
	print STDERR "Read epochs: E=$E\n";
    }    
    if ( /^NUM_TERMS\s/ ) {
	chomp();
	$W = $_;
	$W =~ s/^NUM_TERMS\s+//;
	print STDERR "Read words: W=$W\n";
    }
}
close(I);

open(M,">$TCASTEM.mu");
for (my $e=0; $e<$E; $e++) {
    for (my $k=0; $k<$K; $k++) {
	print M "$e $k $alpha[$k]\n"
    }
}
print STDERR "Wrote '$TCASTEM.mu'\n";

for (my $e=0; $e<$E; $e++) {
    my $tcaname = sprintf("$TCASTEM.phi%03d", $e);
    open(T,">$tcaname");
    for (my $k=0; $k<$K; $k++) {
	@vec = ();
	my $dtmname = sprintf("$DTMDIR/lda-seq/topic-%03d-var-e-log-prob.dat", $k);
	open(D,"<$dtmname");
	#  skip earlier epochs
	for (my $cnt=$W*$e; $cnt>0; $cnt--) {
	    <D>;
	}
	my $tot = 0.0;
	for (my $w=0; $w<$W; $w++) {
	    $vec[$w] = <D>;
	    chomp($vec[$w]);
	    $tot += $vec[$w] = exp($vec[$w]);
	}
	for (my $w=0; $w<$W; $w++) {
	    $vec[$w] /= $tot;
	    print T "$w $k $vec[$w]\n";
	}
	close(D);
    }
    close(T);
    print STDERR "Wrote '$tcaname'\n";
}

system("perl -pi -e 's/^at = .*/at = 0.0/' $TCASTEM.par")==0 or die "Cannot modify $TCASTEM.par\n";
system("perl -pi -e 's/^bb = .*/bb = 0.0/' $TCASTEM.par")==0 or die "Cannot modify $TCASTEM.par\n";
system("perl -pi -e 's/^bt = .*/bt = $conc/' $TCASTEM.par")==0 or die "Cannot modify $TCASTEM.par\n";
print STDERR "Modified $TCASTEM.par\n";

