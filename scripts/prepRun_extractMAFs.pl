#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX"); # no PREs on chrY pass quality filters

open(OUT, ">extractMAFs.batch");
foreach my $i (@chroms){
	print OUT ". ~/.bashrc; ";
	print OUT "cd data/MAFs/; ";
	print OUT "wget https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/maf/$i\.maf.gz; ";
	print OUT "gunzip $i\.maf.gz; ";
	print OUT "zcat ../CREs.bed.gz | grep -P \"^$i\\t.*pCE\" | cut -f 1-3 >pCEs_$i\.bed; ";
	print OUT "maf_parse -g pCEs_$i\.bed $i\.maf >pCEs_$i\.maf; ";

	print OUT "zcat ../CREs.bed.gz | grep -P \"^$i\\t.*cCRE\" | cut -f 1-3 >cCREs_$i\.bed; ";
	print OUT "maf_parse -g cCREs_$i\.bed $i\.maf >cCREs_$i\.maf";
	if($i eq "chr1" | $i eq "chr2"){ # large chromosomes are split into two for phastBias
		print OUT "; maf_parse -S 150000000 -d 1 -r $i\_ $i\.maf";
	}
	print OUT "\n";
}
close(OUT);

# Create jobs:
# module load dSQ
# dsq --job-file extractMAFs.batch --mem-per-cpu 1G --mail-type FAIL,END
