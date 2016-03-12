#!/usr/bin/perl

# generate a topic diagnostics file:
#      hca -v -v -V -e -C0 -r0 ..... STEM RES 2>&1 | ./trim.pl
# but to figure out what the coloums are you will need to
# read the manual

$n = 0;
while ( <> ) {
	if ( / p=/  ) {
		my $withng = 0;
		chomp();
		if ( / ng=/ ) {
			$withng = 1;
 		}
		s/.* p=//;
		s/ ngl=.*//;
		s/ [a-z1]+=/ /g;
		s/,/ /;
		s/%//g;
		print "$n ";
		print;
		$n++;
		$words = <>;
		chomp($words);
		$words =~ s/.*words=//;
		print " $words\n";
	}
}

