#!/usr/bin/perl
use strict;
use warnings;
use utf8;
use List::Util qw(max min);
use Getopt::Std;

my $usage = "
extract-promotor.pl
written by Severin Uebbing                                14 May 2020

Extracts asymmetric intervals from GTF feature starts (strand aware). 
Useful to extract promotor regions.

Usage: extract-promotor.pl -c <chromosome file> [-5 INT -3 INT -i] 
	<infile> > <outfile>
  -c STRING		Chromosome lengths file.
  -5 INT		Upstream extension. Default is 5000 bp.
  -3 INT		Downstream extension. Default is 1000 bp.
  -i			Switches feature name to ENSEMBL gene ID. Default is 
  					official gene name.
  <infile>		Gencode annotations GTF.
  <outfile>		Feature BED file.";

getopts('ic:5:3:');
our($opt_c);
our($opt_5);
our($opt_3);
our($opt_i);
if(@ARGV < 1 | !(defined $opt_c)){die "Error: Input file(s) missing.\n$usage\n";}
my $gencode = $ARGV[0]; # gencode annotation GTF
if(!defined $opt_5){
	$opt_5 = 5000;
}
if(!defined $opt_3){
	$opt_3 = 1000;
}

my %len = ();
open(IN, "<$opt_c");
while(<IN>){
	chomp;
	my @s = split;
	$len{$s[0]} = $s[1];
}
close(IN);

open(IN, "<$gencode");
while(<IN>){
	next if(/^#/);
	my @s = split;
	next if($s[2] !~ /^transcript$/);
	if($s[6] =~ /\+/){
		my $start = max 0, $s[3] - $opt_5 -1;
		my $end = min $len{$s[0]}, $s[3] +$opt_3;
		print "$s[0]\t$start\t$end\t";
		if(!defined $opt_i){
			my @a = split /"/, $s[15];
			print "$a[1]";
		}else{
			my @a = split /"|\./, $s[9];
			print "$a[1]";
		}
	}elsif($s[6] =~ /-/){
		my $end = max 1, $s[4] - $opt_3 -1;
		my $start = min $len{$s[0]}, $s[4] + $opt_5;
		print "$s[0]\t$end\t$start\t";
		if(!defined $opt_i){
			my @a = split /"/, $s[15];
			print "$a[1]";
		}else{
			my @a = split /"|\./, $s[9];
			print "$a[1]";
	}}
	print "\t.\t$s[6]\n";
}
close(IN);
