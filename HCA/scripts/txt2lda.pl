#!/usr/bin/perl -w

#  input a txtbag format and out an ldac format
#     ./txt2lda.pl < STEM.txtbag > STEM.ldac

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use POSIX;
use Getopt::Long;
use Pod::Usage;

my $vocab = "";

my %words = ();
       
GetOptions(
    't|tokens=s'   => \$vocab,
    'man'      => sub {pod2usage(-exitstatus => 0, -verbose => 2)},
    'h|help'       => sub {pod2usage(1)},
);

pod2usage(-message => "ERROR: need input file and stem")
      if ( $#ARGV != 1 );

my $file = shift();
my $ldacfile = shift();

if ( $vocab ne "" ) {
    #  read the map
    open(T,"<$vocab") or die "Cannot open 'vocab': $!\n";
    while ( ($_=<T>) ) {
	chomp($_);
	foreach my $a ( split(/\s+/,$_) ) {
	    $words{$a+0} = 1;
	}
    }
    print STDERR "\n";
    close(T);
    #  create the LDAC tokens file
    open(T,"<$file.tokens") or die "Cannot open '$file.tokens': $!\n";
    open(S,">$ldacfile.tokens") or die "Cannot open 'ldac$file.tokens': $!\n";
    my $line = 0;
    while ( ($_=<T>) ) {
	if ( $words{$line} ) { print S $_; }
	$line++;
    }
    close(T);
    close(S);
}

open(F,"<$file.txtbag") or die "Cannot open '$file.txtbag': $!\n";
open(O,">$ldacfile.ldac") or die "Cannot open '$ldacfile.ldac': $!\n";

<F>;
<F>;
while ( ($_=<F>) ) {
	chomp($_);
	my @a = split(/ /,$_);
	my $len = shift(@a);
	$len = 0;
	my @b = ();
	while ( @a ) {
	  my $i = shift(@a);
	  my $n = shift(@a);
	  if ( !defined($n) ) {
	      exit(0);
	  }
	  if ( $vocab eq "" || $words{$i+0} ) {
	      push(@b,"$i:$n");
	      $len++;
	  }
        }
	print O "$len ";
	print O join(" ",@b);
	print O "\n";
}
close(F);
close(O);

__END__

=head1 NAME
  
txt2lda.pl - convert text bag format to LDA format, optionally thinning tokens

=head1 SYNOPSIS
    
txt2lda.pl [-t tokens] SRCSTEM DESTSTEM

Options:

    SRCSTEM             stem for .txtbag file and optionally .tokens
    DESTSTEM            for outputs
    -h, --help          display help message and exit.
     --man              print man page and exit.
    -t tokens           token file giving tokens to keep

=head1 DESCRIPTION

Convert ".txtbag" format into LdaC format.
Optionally, map a subset of the tokens.  The set to map are
those listed in
I<tokens> file 
as a space separated list of integers.
It is assumed token ids are offset at 0 and numbered sequentially.

=head1 AUTHOR

Wray Buntine

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 Wray Buntine

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.

=cut
