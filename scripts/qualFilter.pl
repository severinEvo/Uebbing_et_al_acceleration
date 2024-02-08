#!/usr/bin/perl
use strict;
use warnings;
use utf8;
use POSIX;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $pos = $infile;
$pos =~ s/.*(chr[XY0-9]+):([0-9]+)-([0-9]+)\.maf/$1\t$2\t$3/;

my %hash = ();
open(IN, $infile);
while(<IN>){
	next if($_ !~ /^s/);
	my @s = split /\s+/;
	my @d = split /\./, $s[1];
	$hash{$d[0]} += $s[3];
}
close(IN);

my $short = POSIX::ceil($hash{'hg38'} * .8); # >= 80% of the human alignment length
my $long = POSIX::floor($hash{'hg38'} * 1.2); # <= 120% of the human alignment length
open(OUT, '>>', $outfile);
foreach my $key (keys %hash){
	if($hash{$key} >= $short & $hash{$key} <= $long){
		print OUT "$pos\t$key\n";
}}
close(OUT);
