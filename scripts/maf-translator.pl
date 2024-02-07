#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my $maf = $ARGV[0]; # PRE MAF file
my $outfile = $maf;
$outfile =~ s/\.maf/.bed/;
open(IN, "<$maf");
while(<IN>){
	next if($_ !~ /^s /);
	my @s = split /\s+/;
	my $start = $s[2];
	my $end = $s[2]+$s[3];
	my @d = split /\./, $s[1];
	my $out = $d[0] . "_" . $outfile;
	open(OUT, ">>$out");
	print OUT "$d[1]\t$start\t$end\t$d[0]\n";
	close(OUT);
}
close(IN);

qx/for i in \$(ls *_$outfile); do sort -k1,1 -k2,2n \$i | bedtools merge -i - -d 1 -c 4 -o distinct >> $outfile; done/;
qx/rm *_$outfile/;

my $hlen;
my %count = ();
open(IN, "<$outfile");
while(<IN>){
	my @s = split /\s+/;
	$count{$s[3]} ++;
	if($_ =~ /hg38/){
		$hlen = $s[2]-$s[1];
}}
close(IN);

my $min = .8 * $hlen;
my $max = 1.2 * $hlen;
my $transfile = $outfile;
$transfile =~ s/\.bed/_trans.bed/;
open(IN, "<$outfile");
open(OUT, ">$transfile");
while(<IN>){
	chomp;
	my @s = split /\s+/;
	my $len = $s[2]-$s[1];
	if($len >$min && $len <$max && $count{$s[3]} == 1){
		print OUT "$_\n";
}}
close(OUT);
