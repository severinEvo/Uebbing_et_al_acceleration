#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX");

open(OUT, ">PRE-MAFs-job.txt");
foreach my $i (@chroms){
	print OUT ". ~/.bashrc; mkdir -p data/MAFs/tmp_mafs/$i ; cd data/MAFs/; ";

	print OUT "LIST=\$(perl -lane 'print \"\$F[0]:\$F[1]-\$F[2]\";' pCEs_$i\.bed); for j in \$LIST; ";
	print OUT "do echo \"\$j\" | perl -lne '\$_ =~ s/(chr[XY0-9]+):([0-9]+)-([0-9]+)/\$1\\t\$2\\t\$3/; print;' >tmp_pCE_$i\.bed; ";
	print OUT "maf_parse -g tmp_pCE_$i\.bed pCEs_$i\.maf >tmp_mafs/$i\/\"\$j\"\\.maf; ";
	print OUT "rm tmp_pCE_$i\.bed; ";
	print OUT "done; ";

	print OUT "LIST=\$(perl -lane 'print \"\$F[0]:\$F[1]-\$F[2]\";' cCREs_$i\.bed); for j in \$LIST; ";
	print OUT "do echo \"\$j\" perl -lne '\$_ =~ s/(chr[XY0-9]+):([0-9]+)-([0-9]+)/\$1\\t\$2\\t\$3/; print;' >tmp_cCRE_$i\.bed; ";
	print OUT "maf_parse -g tmp_cCRE_$i\.bed cCREs_$i\.maf >tmp_mafs/$i\/\"\$j\"\\.maf; ";
	print OUT "rm tmp_cCRE_$i\.bed; ";
	print OUT "done\n";
}
close(OUT);

# Create jobs:
# module load dSQ
# dsq --job-file PRE-MAFs-job.txt --mem-per-cpu 6G --mail-type FAIL,END
