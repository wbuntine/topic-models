#!/usr/bin/perl

#  input a txtbag format and out an ldac format
#     ./txt2lda.pl < STEM.txtbag > STEM.ldac

<>;
<>;
while ( <> ) {
	chomp($_);
	@a = split(/ /,$_);
	$len = shift(@a);
	print "$len";
	while ( @a ) {
	  $i = shift(@a);
	  $n = shift(@a);
	  print " $i:$n";
        }
	print "\n";
}
