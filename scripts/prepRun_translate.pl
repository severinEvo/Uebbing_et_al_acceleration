#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX");

open(OUT, ">PRE-MAFs.batch");
foreach my $i (@chroms){
	print OUT ". ~/.bashrc; mkdir -p data/MAFs/tmp_MAFs/$i ; cd data/MAFs/tmp_MAFs/$i ; ";
	print OUT "module load BEDTools; ";
	print OUT "for i in \$(ls chr1*.maf); do perl ~/whales/200120_120mammals/220920_meme/maf-translator.pl $i ; ";
	print OUT "done\n";
}
close(OUT);

# Create jobs:
# module load dSQ
# dsq --jobfile perl-translate-job.txt --mem-per-cpu 6G --mail-type FAIL,END
