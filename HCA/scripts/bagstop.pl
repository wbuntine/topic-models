#!/usr/bin/perl

#  remove the first #3 tokens/words from the vocab 
#  from dataset #1 to produce dataset #2
#  changing:   .words, .tokens, .txtbag, .colls

$STEMIN = shift();
$STEMOUT = shift();
$chop = shift();

if ( ! -e "$STEMIN.txtbag" ) {
	print STDERR "Intended to modify '$STEMIN.txtbag' and related files\n";
	print STDERR "Couldn't find '$STEMIN.txtbag'\n";
	exit(1);
}

if ( -e "$STEMIN.tokens" ) {
    open(F,"<$STEMIN.tokens") or die
	"Cannot open $STEMIN.tokens";
    open(G,">$STEMOUT.tokens");
    for (my $k=0; $k<$chop; $k++) {
	<F>;
    }
    while ( ($_=<F>) ) {
	print G $_;
    }
    close(F);
    close(G);
}

if ( -e "$STEMIN.words" ) {
    open(F,"<$STEMIN.words") or die
	"Cannot open $STEMIN.words";
    open(G,">$STEMOUT.words");
    for (my $k=0; $k<$chop; $k++) {
	<F>;
    }
    $line = 0;
    while ( ($_=<F>) ) {
	s/^[0-9]+ /$line /;
	print G $_;
	$line++;
    }
    close(F);
    close(G);
}

if ( -e "$STEMIN.txtbag" ) {
    open(F,"<$STEMIN.txtbag") or die
	"Cannot open $STEMIN.txtbag";
    open(G,">$STEMOUT.txtbag");
    $_ = <F>;
    print G $_;
    $_ = <F>;
    chomp();
    $_ -= $chop;
    print G "$_\n";
    while ( ($_=<F>) ) {
	chomp();
	@a = split(/ /, $_);
	my $els = shift(@a);
	for ($k=0; $k<=$#a; $k++) {
	    $a[$k] -= $chop;
	}
	print G "$els " . join(" ",@a) . "\n";
    }
    close(F);
    close(G);
}

if ( -e "$STEMIN.colls" ) {
    open(F,"<$STEMIN.colls") or die
	"Cannot open $STEMIN.colls";
    open(G,">$STEMOUT.colls");
    while ( ($_=<F>) ) {
	chomp();
	@a = split(/ /, $_);
	my $els = shift(@a);
	for ($k=0; $k<=$#a; $k+=3) {
	    $a[$k+2] -= $chop;
	}
	print G "$els " . join(" ",@a) . "\n";
    }
    close(F);
    close(G);
}
