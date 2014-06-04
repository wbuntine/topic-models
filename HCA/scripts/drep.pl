#!/usr/bin/perl

use strict;
use utf8;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use HTML::Entities;

# encoding pragmas follow any includes like "use"
use open ':utf8';
binmode STDIN, ":utf8";
binmode STDERR, ":utf8";

#  constants
my $WORDSPERTOPIC = 20;
my $MINPROB = 0.001;


# globals
my $stopfile;
my %stops = ();
my $verbose = 1;
my $maxN = 10000000;
my @token = ();      #  token/word map, starts at 0
my @docid = ();      #  doc ids, numeric, from .docs file
my @topicxml = ();   #  fragment for each topic
my @docxml = ();     #  fragment for each doc

# dimensions and hyperparameters read from .par file
my $T = 0;
my $W = 0;
my $a0 = 0;
my $b0 = 2;
my $apar = 0;
my $bpar = 0;
my @probs = ();
my $alpha = 1;
my $beta = 1;
my $D = 0;
my $DT = 0;
my $TEST = 0;
my $N = 0;
my @Tt = ();
my @Wt = ();
my @Nd = ();
my $datastem = "";
my $xmlfile = "";

GetOptions(
     'man'       => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
      'stopfile=s' => \$stopfile,
      'xml=s' => \$xmlfile,
      'maxN=i' => \$maxN,
      'words=i' => \$WORDSPERTOPIC,
      'minP=s' => \$MINPROB,
      'v|verbose' => \$verbose,
      'h|help'       => sub {pod2usage(1)}
);

$MINPROB = (0.0+$MINPROB);

pod2usage(-message => "ERROR: need stem")
      if ( $#ARGV != 0 );

my $stem = shift();

# load stopfile
if ( $stopfile ) {
    open(S,"<$stopfile");
    while ( ($_=<S>) ) {
        chomp();
        $stops{$_} = 1;
    }
    close(S);
    if ( $verbose ) {
        print STDERR "Stopwords from '$stopfile' loaded\n";
    }
}

# load parameters from ".par" file
open(T,"<$stem.par");
while( ($_=<T>) ) {
    if ( /^([a-zA-Z0-9]+)\s*= (.*)$/ ) {
	if ( $1 eq "N" ) {
	    $N = $2;
	} elsif ( $1 eq "TEST" ) {
	    $TEST = $2;
	} elsif ( $1 eq "D" ) {
	    $D = $2;
	} elsif ( $1 eq "W" ) {
	    $W = $2;
	} elsif ( $1 eq "T" ) {
	    $T = $2;
	} elsif ( $1 eq "apar" ) {
	    $apar = $2;
	} elsif ( $1 eq "bpar" ) {
	    $bpar = $2;
	} elsif ( $1 eq "a0" ) {
	    $a0 = $2;
	} elsif ( $1 eq "b0" ) {
	    $b0 = $2;
	} elsif ( $1 eq "alpha" ) {
	    $alpha = $2;
	} elsif ( $1 eq "beta" ) {
	    $beta = $2;
	} elsif ( $1 eq "stem" ) {
	    $datastem = $2;
	} elsif ( $1 eq "probs" ) {
	    my $p = $2;
	    $p =~ s/-/0/g;
	    $p =~ s/^\s+//;
	    @probs = split(/ /,$p);
	}
    }
}
close(T);
$DT = $D - $TEST;

if ( $datastem eq "" ) {
    print STDERR "Couldn't find stem in '$stem.par'\n";
    exit(0);
}
if ( $xmlfile eq "" ) {
    $xmlfile = $datastem . ".xml";
}

# load tokens
open(T,"<$datastem.tokens");
my $tc = 0;
while ( ($_=<T>) ) {
    chomp();
    $token[$tc++] = $_;
}
close(T);
if ( $verbose ) {
    print STDERR "Loaded $tc tokens from '$stem.tokens'\n";
}

# load doc info
open(T,"<$datastem.docs");
my $tc = 0;
while ( ($_=<T>) ) {
    chomp();
    my @a = split();
    if ($tc != $a[0] ) {
	print STDERR "Numbers out for docs at $_\n";
	exit(0);
    }
    $docid[$tc] = int($a[2]);
    #print STDERR " $tc/$docid[$tc]";
    $tc++;
}
close(T);
if ( $verbose ) {
    print STDERR "Loaded $tc tokens from '$stem.tokens'\n";
}

open(D,"<$datastem.dit");
while ( ($_=<D>) ) {
    chomp();
    $Nd[$_-1]++;
}
close(D);

#  initialise the topic XMLs
my @txc = ();
for ( my $t=0; $t<$T; $t++) {
    $topicxml[$t] = "<Topic id=\"$t\" prob=\"$probs[$t]\">\n";
    $txc[$t] = 0;
}

# count topic prevalencies and update topic XML
open(T,"sort -k 3 -n -r -t ' ' $stem.wp |");
my $Ncount = 0;
while( ($_=<T>) ) {
    my @a = split();
    if ( $#a==2 ) {
	my $t = $a[1]-1;
	if ( $txc[$t]<$WORDSPERTOPIC ) {
	    if ( ! $stops{$token[$a[0]-1]} ) {
		$topicxml[$t] .= 
		    "  <Word weight=\"$a[2]\">" .
			HTML::Entities::encode_entities($token[$a[0]-1]) . "</Word>\n";
		$txc[$t]++;
	    }
	}
	$Ncount += $a[2];
	$Tt[$t] += $a[2];
	$Wt[$a[0]-1] += $a[2];
    }
}
close(T);
if ( $verbose ) {
    print STDERR "Read $Ncount words from '$stem.wp' out of $N in dataset\n";
}

#  sort topics
my @topicsort = sort { $probs[$b] <=> $probs[$a] } (0..($T-1));

print "<Data id=\"root\">\n";
print "<Documents count=\"$DT\"/>\n";

#  print topics
print "<Topics>\n";
for ( my $st=0; $st<$T; $st++) {
    my $t = $topicsort[$st];
    if ( $txc[$t]>0 && $probs[$t]>$MINPROB ) {
	$topicxml[$t] .= "</Topic>\n";
    } else {
	print STDERR "Topic dropped:  $t  $txc[$t]  p=$probs[$t]\n";
	$topicxml[$t] = "";
    }
    print $topicxml[$t];
}
print "</Topics>\n";

# print documents
open(Z,"<$stem.z");
open(W,"<$datastem.wit");
open(T,"<$xmlfile");
$tc = 0;
while ( defined($_=<T>) && $tc<$maxN ) {
    if ( /^<Document id=\"([0-9]+)\">/ ) {
	if ( $1 == $docid[$tc] ) {
	    my @twords = ();
	    my @tcount = ();
	    print;
	    while ( ($_=<T>) && $_ ne "<Topics/>\n" ) {
		print;
	    }
	    if ( $_ ne "<Topics/>\n" ) {
		print STDERR "Couldn't find Topics element for doc $docid[$tc]\n";
		exit(0);
	    }
	    
	    for (my $z=0; $z<$Nd[$tc]; $z++) {
		my $zin = <Z>-1;
		my $win = <W>-1;
		$tcount[$zin] ++;
		$twords[$zin] .= "  <W>$token[$win]</W>\n";
	    }
	    print " <Topics>\n";
	    for (my $t=0; $t<$T; $t++ ) {
		if ( !defined($tcount[$t]) ) {
		    print "  <Topic id=\"$t\" probability=\"0.0\"/>\n";
		} else {
		    my $pp = ($tcount[$t]+0.000)/$Nd[$tc];
		    print "  <Topic id=\"$t\" probability=\"$pp\">\n";
		    print $twords[$t];
		    print "  </Topic>\n";
		}
	    }
	    print " </Topics>\n";
	    while ( ($_=<T>) && $_ ne "</Document>\n" ) {
		print;
	    }
	    print " </Document>\n";
	    $tc++;
	} else {
	    if ( $verbose ) {
		print STDERR "Skipping doc $tc=$docid[$tc] since got $1, abandoning\n";
		print "</Data>\n";
		close(W);
		close(Z);
		close(T);
		if ( $verbose ) {
		    print STDERR "Loaded $tc xml entries from '$xmlfile'\n";
		}
		exit(0);
	    }
	}
    }
}
close(W);
close(Z);
close(T);
if ( $verbose ) {
    print STDERR "Loaded $tc xml entries from '$xmlfile'\n";
}
print "</Data>\n";

