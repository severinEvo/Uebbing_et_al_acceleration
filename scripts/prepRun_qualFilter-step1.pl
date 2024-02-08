#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX");

open(OUT, ">qualFilter-step1_jobfile.txt");
foreach my $i (@chroms){
	print OUT ". ~/.bashrc; ";
	print OUT "cd data/MAFs/; ";

	print OUT "LIST=\$(perl -lane 'print \"\$F[0]:\$F[1]-\$F[2]\";' cCREs_new_$i\.bed); ";
	print OUT "for j in \$LIST; ";
	print OUT "do perl ../../scripts/qualFilter.pl tmp_mafs/\"\$j\"\\.maf ../qualFilter/cCREs_$i\_use.bed; ";
	print OUT "done; ";

	print OUT "LIST=\$(perl -lane 'print \"\$F[0]:\$F[1]-\$F[2]\";' pCEs_new_$i\.bed); ";
	print OUT "for j in \$LIST; ";
	print OUT "do perl ../../scripts/qualFilter.pl tmp_mafs/\"\$j\"\\.maf ../qualFilter/pCEs_$i\_use.bed; ";
	print OUT "done\n";
}
close(OUT);

# dsq --job-file qualFilter-step1_jobfile.txt --mem-per-cpu 100 --mail-type FAIL,END
