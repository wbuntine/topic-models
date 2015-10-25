#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use utf8;
use POSIX;
use Getopt::Long;
use Pod::Usage;

#  command line args
my $STEM;
my $TRIAL;
#  parameters set somehow
my $VERSION = "hca";
my $TRAIN = 0;
my $TEST = 0;
my $USESVM = 0;
my $D = 0;
my $DIAGNOSTICS = "";
my $GIBBSTRAIN=300;
my $GIBBSEVAL=50;
my $GIBBSEVALBURN=10;
my $PMI = "";
my $BUILDPROB = "";
my $FILE = "";
my $doclass = 0;
my $X = "";
my $Q = "";
my $usephi = 0;
my $HOLD = "-hdoc,4";

my @tag = ();
my %command = ();

#  echo command before running and print error on death
sub mysystem() {
      my $comm = shift();
      my $type = shift();
      if ( !defined($type) ) { $type = ""; }
      my $file = $TRIAL;
      print STDERR "Running: $comm\n";
      if ( $comm =~ / ($TRIAL[A_Z])[^[:alpha:]]/ ) {
	  $file = $1;
      }
      system($comm)==0
	  or die "Cannot run $type on $file: " . substr($comm,0,15);
}

#  pick up topic words+details for this call and output to ".tpk"
sub picklast() {
    my $hcamatch = shift();
    my $stem = shift();
    local *L;
    open(L, "<$stem.log");
    while ( ($_=<L>) ) {
	if ( /hca .*$hcamatch/ ) {
	    last;
	}
    }
    if ( !defined($_) ) {
	print STDERR "No topic detail in '$stem.log' matching '$hcamatch'\n";
	close(L);
	return;
    }
    while ( ($_=<L>) ) {
	if ( /^Topic [0-9]+\/0 p=/ ) {
	    last;
	}
    }
    if ( !defined($_) ) {
	print STDERR "No topic detail in '$stem.log' matching '$hcamatch'\n";
	close(L);
	return;
    }
    local *T;
    open(T, ">$stem.tpk");
    print T;
    while ( ($_=<L>) ) {
	if ( /^probs = / ) {
	    last;
	}
	print T;
    }
    close(L);
    close(T);
}

#  grab para. from .srcpar file
sub getpar() {
      my $stem = shift();
      my $par = shift();
      open(SRCPAR, "grep $par= $stem.srcpar |");
      my $fc = <SRCPAR>;
      close(SRCPAR);
      if ( !defined($fc) ) {
        print STDERR "Cannot find par '$par' in '$stem.srcpar'\n";
        exit(1);
      }
      chomp($fc);
      $fc =~ s/^.*=//;
      if ( $fc eq "" ) {
        print STDERR "Cannot find valid par '$par' in '$stem.srcpar'\n";
        exit(1);
      }
      return $fc;
}

sub dotag {
    my @vars = @_;
    print "The vars for tag are: @vars\n";
    if ( $vars[1] =~ /^([a-zA-Z\-0-9]+),(.+)$/ ) {
	$command{$1} = $2;
    } else {
	pod2usage(-message => "ERROR: bad argument to -tag")
    }
}

# print STDERR "COMM: $0 ". (join " ", grep { if ( / / ) { $_ = "\'$_\'"; } else { $_; } } @ARGV) . "\n";

GetOptions(
    'gibbstrain=i'     => \$GIBBSTRAIN,
    'gibbseval=i'     => \$GIBBSEVAL,
    'gibbsburn=i'     => \$GIBBSEVALBURN,
    'version=s'     => \$VERSION,
    'tag=s'     => \&dotag,
    'hold=s'     => sub { $_=shift(); $HOLD="-h$_"; },
    'phi!'     => \$usephi,
    'diagnostics=s' => \$DIAGNOSTICS,
    'X=s'     => \$X,
    'q=i'     => \$Q,
    'A=s'     => \$command{"A"},
    'B=s'     => \$command{"B"},
    'C=s'     => \$command{"C"},
    'D=s'     => \$command{"D"},
    'E=s'     => \$command{"E"},
    'F=s'     => \$command{"F"},
    'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
    'help'   => sub {pod2usage(1)},
    ) or exit(1);

pod2usage(-message => "ERROR: need command line: STEM TRIAL")
      if ( $#ARGV != 1 );

$STEM=shift();
$TRIAL=shift();

#  strip away null commands
foreach my $k ( keys %command ) {
	if ( $command{$k} eq "" ) {
		delete $command{$k};
	}
}
@tag = sort keys %command;

if ( $TRIAL eq $STEM ) {
	print STDERR "Trial name '$TRIAL' should be different from stem '$STEM'\n";
	exit(1);
}

if ( -e "$STEM.pmi" || -e "$STEM.pmi.gz") {
	$PMI = "-p ";
}
if ( $Q ne "" ) {
    $Q = "-q " . $Q;
}
if ( $USESVM && -e "$STEM.class" ) {
	$doclass = 1;
	$BUILDPROB = "-X -ltestprob,1,2 ";
}
if ( -e "$STEM.ldac" ) {
        $FILE = "-f ldac "
} elsif ( -e "$STEM.wit" ) {
        $FILE = "-f witdit "
} elsif ( -e "$STEM.txtbag" ) {
        $FILE = "-f bag "
}

$TEST = &getpar($STEM,"testdocs");
$D = &getpar($STEM,"documents");
if ( $D<=$TEST ) {
	print STDERR "Test count=$TEST, doc count=$D!!\n";
	exit(1);
}
$TRAIN=$D-$TEST;

print STDERR "Running tests:   STEM=$STEM TRAIL=$TRIAL + (" . join(", ",@tag) . ") TEST=$TEST $X\n";
foreach my $A ( sort @tag ) {
	print STDERR "  trial $A :: $X $command{$A} \n";
}
#  clean up
foreach my $A ( @tag ) {
	unlink glob "${TRIAL}$A.*";
}

#   build models
#   but switch off likelihood calcs at end
foreach my $A ( @tag ) {
  &mysystem("$VERSION $Q $X $FILE -t$TRAIN  -Llike,0,0 -v $command{$A} -C$GIBBSTRAIN $STEM ${TRIAL}$A", "hca train");
}
#   do doc-topic probability estimation over $GIBBSEVAL samples
#   and switch off likelihood calcs at end
foreach my $A ( @tag ) {
  &mysystem("$VERSION $Q $FILE -C$GIBBSEVAL $BUILDPROB -Llike,0,0 -lphi,1,2 $DIAGNOSTICS -v -v -r 0 -T$TEST $STEM ${TRIAL}$A","hca eval");
}
if ( $doclass ) {
     #   do holdout perp calc over $GIBBSEVAL cycles including a $GIBBSEVALBURN burnin
     #   compute PMI if a .pmi file exists
     #   generate SVM train and test files and run SVM (NB. .testprob files generated by previous runs using $BUILDPROB)
   foreach my $A ( @tag ) {
     &mysystem("$VERSION $FILE $PMI -C0 -X -v $HOLD -Llike,$GIBBSEVAL,$GIBBSEVALBURN -v -V -r 0 -T$TEST $STEM ${TRIAL}$A","hca test");
     if ( $USESVM ) {
	if ( "${TRIAL}$A.testprob" ) {
           system("head -$TRAIN ${TRIAL}$A.testprob | sed -e 's/^[0-9]*: //' > ${TRIAL}$A.lstrain");
           system("tail -n $TEST ${TRIAL}$A.testprob | sed -e 's/^[0-9]*: //' > ${TRIAL}$A.lstest");
           &mysystem("svm-train ${TRIAL}$A.lstrain ${TRIAL}$A.lsmodel > ${TRIAL}$A.lslog 2>&1");
           &mysystem("svm-predict ${TRIAL}$A.lstrain ${TRIAL}$A.lsmodel ${TRIAL}$A.lsrun >> ${TRIAL}$A.lslog 2>&1");
        } else {
	  print STDERR "Cannot locate SVM data input file '${TRIAL}$A.tprop'\n";
        }
     }
   }
} else {
   foreach my $A ( @tag ) {
     &mysystem("$VERSION $FILE $PMI -C0 -v $HOLD -Llike,$GIBBSEVAL,$GIBBSEVALBURN -v -V -r 0 -T$TEST $STEM ${TRIAL}$A","hca test");
   }
}
if ( $usephi ) {
   foreach my $A ( @tag ) {
     &mysystem("$VERSION $X $FILE -C0 -v $HOLD -Llike,$GIBBSEVAL,$GIBBSEVALBURN -v -r phi -T$TEST $STEM ${TRIAL}$A","hca test");
   }
}
__END__

=head1 NAME

runtest - run sequence of hca tests.

=head1 SYNOPSIS

runtest STEM TRIAL

Options:

    STEM                stem of input data set
    TRIAL               stem form trial output
    -A ARGS             addition arguments for trial A
    -B ARGS             addition arguments for trial B
    -C/-D/.../-F  ARGS       optional 3rd, 4th, ... 6th trials
    --diagnostics S     add string S as args to diagnostics run of hca
    --gibbsburn C        burnin cycles as part of --gibbseval (deflt=10)
    --gibbseval C        major gibbs cycles with 'hca -l'(deflt=50)
    --gibbstrain C       major gibbs cycles for training (deflt=300)
    -h, --help          display help message and exit.
    --tag TAG,COMMAN    optional trial labelled TAG using command COMMAND
                        TAG is any sequence of letters, numbers or dashes.
    --version HCA       use this executable instead of the default 'hca'
    -X ARGS             args given to every training run
    --man               print man page and exit.

=head1 DESCRIPTION

Runs a sequence of tests with various options under trials under tags
A, B, C, D, E or F, or optionally a user defined tag.
The 
F<STEM.srcpar> file
must have an entry "testdocs=" giving the number of test documents.
If 
F<STEM.class> file exists 
then also run libsvm in default mode
on the document topic proportions, as well as doing standard
Bayes classifier estimate of class based on these.  
If a 
F<STEM.pmi> or 
F<STEM.pmi.gz> then record PMI for the topics.
Perplexity is done with document completion.

Outputs/results go to various extensions of
F<$(TRIAL)A> for trial A
and
F<$(TRIAL)B> for trial B, etc.

=head1 EXAMPLE

	runtest -X -K20 -A '-A1 -B1' -B '-B1' -C ' '  data/ch /tmp/c
	runtest -X --tag AA,'-A1 '  data/ch /tmp/c

First one runs with 100 topics and three trials, A,B,C. 
Output will be to files with stem "cA", "cB" and "cC".
Second one does an AA trial with 20 topics
with command options '-A1 '. Output will be to files with stem "cAA".

=head1 SEE ALSO

I<hca>(1).

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 Wray Buntine

This programme is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.
