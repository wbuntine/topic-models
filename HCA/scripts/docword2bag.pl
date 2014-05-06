#!/usr/bin/perl -w

use POSIX;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;

my $LDAC = 0;

GetOptions(
     'man'       => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
     'ldac'      => \$LDAC,
#      'fields=i' => \$MAXREP,
      'h|help'       => sub {pod2usage(1)}
);

pod2usage(-message => "ERROR: need stem")
      if ( $#ARGV != 0 );

my $STEM = shift();

# copy token file

if ( -f "vocab.$STEM.txt" ) {
	system("cp vocab.$STEM.txt $STEM.tokens");
} else {
	die "Cannot find 'vocab.$STEM.txt' for read";
}

# create txtbag file

if ( -f "docword.$STEM.txt" ) {
	open(DW,"<docword.$STEM.txt") or die "Cannot open 'docword.$STEM.txt' for read";
} elsif ( -f "docword.$STEM.txt.gz" ) {
	open(DW,"zcat docword.$STEM.txt.gz |") or die "Cannot open 'docword.$STEM.txt.gz' for read";
} else {
	die "Cannot find 'docword.$STEM.txt' for read";
}
if ( $LDAC == 0 ) {
   open(TB,">$STEM.txtbag");
   print STDERR "Converting docword format to txtbag format\n";
} else {
   open(TB,">$STEM.ldac");
   print STDERR "Converting docword format to ldac format\n";

}
print STDERR "==========================================\n";

my $DOCS = int(<DW>);
my $WORDS = int(<DW>);
if ( $LDAC == 0 ) {
  print TB "$DOCS\n$WORDS\n";
}
<DW>;

$doc = 1;
%cnt = ();
$nl = 0;
$lf = 0;
$lc = 0;

sub num() {  return $a-$b; }

print STDERR "Done documents:";
while ( ($_=<DW>) ) {
  ($d,$f,$c) = split();
  if ( $d==$doc+1 )  {
      if ( $doc>0 && ($doc % 1000)==0 ) {
	  print STDERR "$doc, ";
      }
      $doc++;
      if ( $lc>0 ) {
	  $cnt{$lf} = $lc;
	  $nl ++;
      }
      print TB "$nl ";
      foreach $k ( sort { $a - $b }  (keys %cnt) ) {
          if ( $LDAC == 0 ) {
	  	print TB "$k $cnt{$k} ";
	  } else {
		print TB "$k:$cnt{$k} ";
	  }
      }
      print TB "\n";
      $nl = 0;
      %cnt = ();
      $lf = 0;
      $lc = 0;
  }
  if ( $f<=0 || $f> $WORDS ) {
   die "Bad feature in: $_";
  }
  if ( $c<=0 ) {
    die "Bad doc, got negative count in: $_";
  }
  if ( $d==$doc ) {
    if ( $f==$lf ) {
      $lc += $c;
    } else {
      if ( $lc>0 ) {
	$cnt{$lf} = $lc;
	$nl++;
      }
      $lc = $c;
      $lf = $f;
    }
  } elsif ( $d>$doc && $d<$doc+10 ) {
    #  empty docs
    for (  ; $doc<$d; $doc++ ) {
      print TB "0\n";
    }
  } else {
    die "Bad doc, got $d, expecting $doc in $_";
  }
}
if ( $lc>0 ) {
  $cnt{$lf} = $lc;     
  $nl ++;
}
print TB "$nl ";
foreach $k ( sort(keys %cnt) ) {
  print TB "$k $cnt{$k} ";
}
print TB "\n";
print STDERR "\n";

close(TB);
close(DW);

exit 0;

__END__

=head1 NAME

docword2bag -- convert UCI docword format to a text bag format

=head1 SYNOPSIS

docword2bag STEM

  Options:
    -h, --help          display help message and exit.
    --ldac              output the LDAC format.
    --man               print man page and exit.

=head1 DESCRIPTION

Expects to find vocab.STEM.txt and docword.STEM.txt or docword.STEM.txt.gz
in the UCI data format in the current directory.  
These are converted to the input required for
DCA's ".txtbag" format or the LDAC format.

This a slow Perl script, so conversion is not fast.

=head1 SEE ALSO

I<linkTables>(1),

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut


