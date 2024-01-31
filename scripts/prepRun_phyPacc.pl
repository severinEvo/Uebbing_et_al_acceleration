#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chroms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX");
my @clades = ("pan","hominini","homininae","hominidae","hominoidea","papionini",
	"cercopithecinae","presbytini","colobini","colobinae","cercopithecidae","catarrhini",
	"callitrichoidea","cebidae","platyrrhini","simiiformes","haplorrhini","lemuroidea",
	"strepsirrhini","primates","primatomorpha","euarchonta","arvicolinae-cricetinae",
	"cricetidae","mus","murinae","muridae","eumuroida","muroidea","myodonta","castorimorpha",
	"phiomorpha","octodontoidea-chinchilloidea","caviomorpha","hystricomorpha",
	"myodonta-castorimorpha-hystricomorpha","sciuromorpha","rodentia","lagomorpha","glires",
	"euarchontoglires","camelus","camelidae","delphinidae","delphinoidea",
	"delphinoidea-inoidea","odontoceti","mysticeti","cetacea","bovina","bovini","ovis",
	"caprinae","caprini","bovidae","cervidae","pecora","cetruminantia","artiofabula",
	"cetartiodactyla","equus","equidae","perissodactyla","felinae","panthera","felidae",
	"canidae","mustelidae","musteloidea","ursidae","phocidae","pinnipedia","arctoidea",
	"caniformia","carnivora","pholidota","ferae","pteropus","pteropodidae","yinochiroptera",
	"yinpterochiroptera","vespertilionoidea","yangochiroptera","chiroptera","scrotifera",
	"erinaceidae-soricidae","eulipotyphla","laurasiatheria","boreoeutheria","tethytheria",
	"paenungulata","afrosoricida","afroinsectivora","afroinsectiphilia","afrotheria",
	"xenarthra","atlantogenata","eutheria","australidelphia","marsupialia");

foreach my $i (@chroms){
	open(OUT, ">phyPacc_$i\_jobfile.batch");
	foreach my $j (@clades){
		print OUT ". ~/.bashrc; ";
		print OUT "cd data/phyloP_out/; ";
		print OUT "phyloP -i MAF -m LRT -f ../MAFs/cCREs_$i\.bed -o ACC -B $j ../../annotations/tree_hg38_120mammals_named.mod ../MAFs/cCREs_$i\.maf >$i\_$j\_cCREs_phyP_ACC.out\n";
		print OUT ". ~/.bashrc; ";
		print OUT "cd data/phyloP_out/; ";
		print OUT "phyloP -i MAF -m LRT -f ../MAFs/pCEs_$i\.bed -o ACC -B $j ../../annotations/tree_hg38_120mammals_named.mod ../MAFs/pCEs_$i\.maf >$i\_$j\_pCEs_phyP_ACC.out\n";
	}
	close(OUT);
}

# Create jobs:
# module load dSQ
# dsq --job-file phyPacc_chr1_jobfile.batch --mem-per-cpu 6G --mail-type FAIL,END
# dsq --job-file phyPacc_chr2_jobfile.batch --mem-per-cpu 6G --mail-type FAIL,END
# dsq --job-file phyPacc_chr3_jobfile.batch --mem-per-cpu 4G --mail-type FAIL,END
# dsq --job-file phyPacc_chr4_jobfile.batch --mem-per-cpu 4G --mail-type FAIL,END
# dsq --job-file phyPacc_chr5_jobfile.batch --mem-per-cpu 4G --mail-type FAIL,END
# dsq --job-file phyPacc_chr6_jobfile.batch --mem-per-cpu 4G --mail-type FAIL,END
# dsq --job-file phyPacc_chr7_jobfile.batch --mem-per-cpu 4G --mail-type FAIL,END
# dsq --job-file phyPacc_chr8_jobfile.batch --mem-per-cpu 3G --mail-type FAIL,END
# dsq --job-file phyPacc_chr9_jobfile.batch --mem-per-cpu 3G --mail-type FAIL,END
# dsq --job-file phyPacc_chr10_jobfile.batch --mem-per-cpu 3G --mail-type FAIL,END
# dsq --job-file phyPacc_chr11_jobfile.batch --mem-per-cpu 3G --mail-type FAIL,END
# dsq --job-file phyPacc_chr12_jobfile.batch --mem-per-cpu 3G --mail-type FAIL,END
# dsq --job-file phyPacc_chr13_jobfile.batch --mem-per-cpu 2G --mail-type FAIL,END
# dsq --job-file phyPacc_chr14_jobfile.batch --mem-per-cpu 2G --mail-type FAIL,END
# dsq --job-file phyPacc_chr15_jobfile.batch --mem-per-cpu 2G --mail-type FAIL,END
# dsq --job-file phyPacc_chr16_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chr17_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chr18_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chr19_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chr20_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chr21_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chr22_jobfile.batch --mem-per-cpu 1G --mail-type FAIL,END
# dsq --job-file phyPacc_chrX_jobfile.batch --mem-per-cpu 3G --mail-type FAIL,END
