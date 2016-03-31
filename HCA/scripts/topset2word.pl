#!/usr/bin/perl
#
#   A man entry for this script is given at the end using POD.
#   The script class the following commands:
#        wordcloud_cli.py,  dot,

use Getopt::Long;
use Pod::Usage;
use POSIX;

$TOPICS = 0;

$noimages = 0;

#   hue for background topic
$BGHUE = 200;
#  max. hue for other topics (0,...,MAX)
$MAXHUE = 45;
#  image dims
$IMAGEWIDTH = 400;
$IMAGEHEIGHT = 200;
#  only record topcor values over this
$CCMIN = 0.04;
#  allow ($CCFACT * #topics) arcs 
$CCFACT = 3;
#  used to discretise the sizes for tag cloud
$SIZEFACT = 4;
#  ignore topics small than this times the max
$PROPFACT = 0.03;
$MAXWORDS = 10;
#  max width for arcs
$MAXWIDTH = 12;
$EDGEWEIGHT = 1;
$DOT = "";
$LANG = "svg";

$COMMAND = join(" ",@ARGV);
GetOptions(
    'man'       => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
    'maxwidth=i' => \$MAXWIDTH,
    'ccfact=s' => \$CCFACT,
    'ccmin=s' => \$CCMIN,
    'noimages!' => \$noimages,
    'lang=s' => \$LANG,
    'dot=s' => \$DOT,
    'maxwords=i' => \$MAXWORDS,
    'sizefact=i' => \$SIZEFACT,
    'propfact=s' => \$PROPFACT,
    'edgeweight=s' => \$EDGEWEIGHT,
    'h|help'       => sub {pod2usage(1)}
    );

$STEM = shift();
$RES = shift();

if ( !defined($RES) || $RES =~ /\// ) {
	print STDERR "Result stem second argument, must have no directory\n";
	exit(0);
}

sub refact() {
    my $pp = int(shift()*shift()+2.3);
    return $pp;
}

open(F,"<$STEM.topset");
{
    #  first check we have the right version of the file
    $_ = <F>;
    my @a = split();
    if ( $a[0] ne "#topic" || $a[3] ne "prop" || $a[6] ne "eff-words" ) {
	print STDERR "Topic header line wrong in '$STEM.topset'\n";
	exit(0);
    }
    $_ = <F>;
    my @a = split();
    if ( $a[0] ne "#word" || $a[4] ne "count" ||
	 $a[7] ne "df" ) {
	print STDERR "Word header line wrong in '$STEM.topset'\n";
	exit(0);
    }
}
$minbgprop = 10000;
while ( ($_=<F>) ) {
    if ( /^topic ([\-0-9]+) [\-0-9]+ ([0-9\.]+) [0-9\.]+ [0-9\.]+ ([0-9\.]+)/ ) {
	if ( $1>=0 ) {
	    $efw[$1] = $3;
	    $prop[$1] = $2;
	    if ( $2>$maxprop ) {
		$maxprop = $2;
	    }
	}
    }  
    if ( /^word ([\-0-9]+) / ) {
	chomp();
	$topic = $1;
	@s = split(/ /, $_);
	#  print a minimum of 4 words, 
	#  but otherwise shrink count if a small topic
	if ( $s[3]>=$MAXWORDS ) {
	    # $s[3]>= $MAXWORDS*1.001*sqrt($prop[$topic]/$maxprop) ) {
	    #   assumes $s[3] is rank
	    next;
	}
	$cnt = int($s[4]);	
	$cnt = sqrt($cnt);
	if ( $topic>=0 ) {
	    $df = $s[7];
	    $rank = $cnt / ($df + 0.001);
	    $data[$topic]{$s[9]} = "$cnt,$df,$rank";
	    # print "$s[9] $cnt,$df,$rank\n";
	    if ( $rank>$maxrank[$topic] ) {
		$maxrank[$topic] = $rank;
	    }
	    if ( $cnt>$maxcnt[$topic] ) {
		$maxcnt[$topic] = $cnt;
	    }
	    if ( !defined($minrank[$topic]) || $rank<$minrank[$topic] ) {
		$minrank[$topic] = $rank;
	    }
	} else {
	    $databg{$s[9]} = "$cnt,$s[5]";
	    if ( $s[5]>$maxbgprop ) {
		$maxbgprop = $s[5];
	    }	
	    if ( $s[5]<$minbgprop ) {
		$minbgprop = $s[5];
	    }	
            if ( $cnt<$maxbgcnt ) {
                $maxbgcnt = $cnt;
            }

	}
	if ( $topic>=$TOPICS ) {
	    $TOPICS = $topic+1;
	}
    }
}
close(F);
print STDERR "Got $TOPICS topics\n";

if ( ! $noimages ) {
    my $printtop = 0;
    print STDERR "Create images: ";
    mkdir($RES);
    #  compute scales for efw
    $maxefw = 0;
    $minefw = 100000;
    for (my $t=0; $t<$TOPICS; $t++) {
	if ( $prop[$t]<$maxprop*$PROPFACT ) {
	    next;
	}
	my $efw = $efw[$t];
	if ( $efw>$maxefw ) {
	    $maxefw = $efw;
	}
	if ( $efw<$minefw ) {
	    $minefw = $efw;
	}
    }
    #  now print topics
    for (my $t=0; $t<$TOPICS; $t++) {
	if ( $prop[$t]<$maxprop*$PROPFACT ) {
	    next;
	}
	$printtop++;
	open(F,">$RES.txt");
	# print "$t";
	my $href = $data[$t];
	foreach my $k ( keys(%$href) ) {
	    my ($cnt,$df,$rank) = split(/,/, $data[$t]{$k});
	    $rank = $rank / ($maxrank[$t] + 0.0001);
	    $cnt = $cnt / ($maxcnt[$t] + 0.0001);
	    #  message to wordcloud.py
	    #     word,text-size-factor,text-whiteness-factor
	    print F " $k,$rank,$cnt";
	}
	print F "\n";
	close(F);
	#   compute hue
	my $efw = ($efw[$t]-$minefw) / ($maxefw-$minefw);
	my $hue = int($efw*$MAXHUE);
	#  get image dimensions
	my $pp = &refact($prop[$t]/$maxprop,$SIZEFACT);
	# print "$t $prop[$t] $pp\n";
	my $S1 = int(($IMAGEWIDTH * $pp) / ($SIZEFACT+2));
	my $S2 = int(($IMAGEHEIGHT * $pp) / ($SIZEFACT+2));
	#  run
	system("wordcloud_cli.py --text $RES.txt --imagefile $RES/$t.png --background=hue=$hue --width=$S1 --height=$S2")==0
	    or die "Failure of executable 'wordcloud_cli.py'\n";
	print STDERR " $t";
	unlink("$RES.txt");
    }
    if ( %databg ) {
	#  yes, there is a background topic
	open(F,">$RES.txt");
	foreach my $k ( keys(%databg) ) {
	    my ($cnt,$rank) = split(/,/, $databg{$k});
	    if ( $maxbgcnt == 0 ) {
		$cnt = 1;
	    } else {
		$cnt = $cnt / ($maxbgcnt + 0.0001);
	    }
	    $rank = ($rank - $minbgprop) / $maxbgprop;
	    #  message to wordcloud.py
	    #     word,text-size-factor,text-whiteness-factor
	    print F " $k,$cnt,0";
	}
	print F "\n";
	close(F);
	system("wordcloud_cli.py --text $RES.txt --imagefile $RES/BG.png --background=hue=$BGHUE --width=$IMAGEWIDTH --height=$IMAGEHEIGHT")==0
	    or die "Failure of executable 'wordcloud_cli.py'\n";
	unlink("$RES.txt");
	print STDERR " BG";
    }
    print STDERR "\n";
    print STDERR "Printed $printtop topics\n";
}

print STDERR "Create file\n";
open(D,">$RES.dot");
my $printtop = 0;
my $label = $STEM;
#  remove dubious parts of label
$label =~ s/[-\/\#\.]//g;
print D "// Graph generated using \"$COMMAND\"\n";
print D "digraph $label {\nbgcolor=black\nrankdir=LR\nedge [dir=none]\n" .
    "splines=true\n";
for (my $t=0; $t<$TOPICS; $t++) {
    if ( $prop[$t]<$maxprop*$PROPFACT ) {
	next;
    }
    $printtop++;
    $intop[$t]++;
    print D "n$t [image=\"$RES/$t.png\"]\n";
}
if ( %databg ) {
    print D "n$t [image=\"$RES/BG.png\"]\n";
}

#read topic correlations
my @topcor;
my @cor = ();
my $cors = 0;
my $maxcor = 0;
open(C,"<$STEM.topcor") or die "no $STEM.topcor";
while ( ($_=<C>) ) {
    if ( /^([0-9]+) ([0-9]+) ([0-9\.]+)/ ) {
	if ( $intop[$1] && $intop[$2] && $3>$CCMIN ) {
	    $cors++;
	    $topcor[$1][$2] = $topcor[$2][$1] = $3;
	    push(@cor,$3);
	}
    }
}
close(C);
@cor = sort { $b <=> $a } @cor;
$maxcor = $cor[0];
print STDERR "Read $cors entries from $STEM.topcor, max=$maxcor\n";

if ( $CCFACT * $printtop >= $#cor ) {
    # use all correlations
} else {
    $CCMIN = $cor[int($CCFACT * $printtop)];
}
my $printcor = 0;
for (my $t=0; $t<$TOPICS; $t++) {
    if ( $prop[$t]<$maxprop*$PROPFACT ) {
	next;
    }
    for (my $s=0; $s<$t; $s++) {
	if ( $prop[$s]<$maxprop*$PROPFACT ) {
	    next;
	}
	if ( $topcor[$t][$s]>$CCMIN ) {
	    my $width = &refact($topcor[$t][$s]/$maxcor,$MAXWIDTH);
	    # print STDERR "$t -> $s $topcor[$t][$s] $maxcor $width\n";
	    $printcor++;
	    my $ew = $width*$EDGEWEIGHT;
	    print D "n$t -> n$s [color=white,penwidth=$width,weight=$ew]\n";
	}
    }
}
print D "}\n";
print STDERR "printed $printcor correlations\n";
print STDERR "Running: dot -T$LANG -O $DOT $RES.dot\n";
system("dot -T$LANG -O $DOT $RES.dot")==0 
    or die "Cannot find executable 'dot'\n";
close(D);


__END__

=pod

=head1 NAME

topset2word.pl - Create graphic from hca STEM.topset/.topcor files

=head1 SYNOPSIS

topset2word.pl [options] STEM RES

  Options:
   --ccfact=f          stop adding smaller correlations once (f*#topics) got
   --ccmin=f           ignore correlations less than this
   --dot=s             pass on string to dot
   --edgeweight=f      edges have proportional weight to penwidth, by f
   --help              brief help message
   --lang=s            sends to "-T" option in dot ("svg","pdf","png", ...)
   --man               full documentation
   --maxwidth=w        maximum penwidth for arcs
   --maxwords=i        max words per topic
   --noimage           don't compute images, assume OK
   --propfact=i        ignore topics with proportion this much less than max
   --sizefact=i        tag clouds are in scales of 2, 3, ... sizefact+2

=head1 DESCRIPTION

B<This program> reads details from the 
F<STEM.topset> file 
(diagnostic data on topics and their top words)
and the 
F<STEM.topcor> file 
(topic correlations)
and creates an image using 
L<dot> and 
L<wordcloud> .
Results are placed in 
F<RES.dot.svg>, unless the --lang flag is used.

To create the two input files, run 
I<hca> with the right arguments, such as:

=over 4

hca -v -v -V -V -r0 -C0 DATA STEM

=back

Try different options such as
changing the edge weights with
"--edgeweight", changing the formatter to 
L<dot> using "--dot -Kfdp",
formatter, 
and generate different output files
using "--lang".
Also, change the number of links with "--ccfact".

=head1 REQUIRED

Needs 
I<graphviz> installed to run 
I<dot> command.
Also, uses a specially  modified version of
I<wordcloud.py> that must replace the same in
Andreas Mueller's
I<wordcloud> Python system
to run 
I<wordcloud_cli.py> .

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2016 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut

