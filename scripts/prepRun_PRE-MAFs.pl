#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX");

open(OUT, ">PRE-MAFs-job.txt");
foreach my $i (@chroms){
	print OUT ". ~/.bashrc; mkdir -p data/MAFs/tmp_mafs/$i ; cd data/MAFs/tmp_mafs/$i ; ";

	print OUT "LIST=\$(perl -lane 'print \"\$F[0]:\$F[1]-\$F[2]\";' ../../pCEs_$i\.bed); for j in \$LIST; ";
	print OUT "do echo \"\$j\" | perl -lne '\$_ =~ s/(chr[XY0-9]+):([0-9]+)-([0-9]+)/\$1\\t\$2\\t\$3/; print;' >\"\$j\"\\.bed; ";
	print OUT "maf_parse -g \"\$j\"\\.bed pCEs_$i\.maf >\"\$j\"\\.maf; ";
	print OUT "done; ";

	print OUT "LIST=\$(perl -lane 'print \"\$F[0]:\$F[1]-\$F[2]\";' ../../cCREs_$i\.bed); for j in \$LIST; ";
	print OUT "do echo \"\$j\" | perl -lne '\$_ =~ s/(chr[XY0-9]+):([0-9]+)-([0-9]+)/\$1\\t\$2\\t\$3/; print;' >\"\$j\"\\.bed; ";
	print OUT "maf_parse -g \"\$j\"\\.bed cCREs_$i\.maf >\"\$j\"\\.maf; ";
	print OUT "done\n";
}
close(OUT);

# Create jobs:
# module load dSQ
# dsq --job-file PRE-MAFs-job.txt --mem-per-cpu 6G --mail-type FAIL,END
