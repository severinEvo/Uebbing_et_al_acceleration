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
	open(OUT, ">phastBias_$i\_jobfile.txt");
	foreach my $j (@clades){
		if($i eq "chr1" | $i eq "chr2"){
			print OUT ". ~/.bashrc; ";
			print OUT "cd data/phastBias_out/; ";
			print OUT "phastBias --output-tracts $i\_1_$j\_phastBias.gff ../MAFs/$i\_1.maf ../../annotations/tree_hg38_120mammals_named.mod $j >$i\_1_$j\_phastBias.wig; ";
			print OUT "perl -lne '\$_ =~ s/hg38/$i\/; print;' $i\_1_$j\_phastBias.gff >tmp_1_$i\_$j && mv tmp_1_$i\_$j $i\_1_$j\_phastBias.gff; ";
			print OUT "phastBias --output-tracts $i\_2_$j\_phastBias.gff ../MAFs/$i\_2.maf ../../annotations/tree_hg38_120mammals_named.mod $j >$i\_2_$j\_phastBias.wig; ";
			print OUT "perl -lne '\$_ =~ s/hg38/$i\/; print;' $i\_2_$j\_phastBias.gff >tmp_2_$i\_$j && mv tmp_2_$i\_$j $i\_2_$j\_phastBias.gff; ";
			print OUT "cat $i\_1_$j\_phastBias.gff $i\_2_$j\_phastBias.gff >tmp_$i\_$j\_phastBias.gff; ";
			print OUT "rm $i\_1_$j\_phastBias.wig $i\_2_$j\_phastBias.wig\n";
		}else{
			print OUT ". ~/.bashrc; ";
			print OUT "cd data/phastBias_out/; ";
			print OUT "phastBias --output-tracts $i\_$j\_phastBias.gff ../MAFs/$i\.maf ../../annotations/tree_hg38_120mammals_named.mod $j >$i\_$j\_phastBias.wig; ";
			print OUT "perl -lne '\$_ =~ s/hg38/$i\/; print;' $i\_$j\_phastBias.gff >tmp_$i\_$j && mv tmp_$i\_$j $i\_$j\_phastBias.gff; ";
			print OUT "rm $i\_$j\_phastBias.wig\n";
	}}
	close(OUT);
}

# Create jobs:
# module load dSQ
# dsq --job-file phastBias_chr1_jobfile.txt --mem-per-cpu 90G --mail-type FAIL,END
# dsq --job-file phastBias_chr2_jobfile.txt --mem-per-cpu 100G --mail-type FAIL,END
# dsq --job-file phastBias_chr3_jobfile.txt --mem-per-cpu 90G --mail-type FAIL,END
# dsq --job-file phastBias_chr4_jobfile.txt --mem-per-cpu 70G --mail-type FAIL,END
# dsq --job-file phastBias_chr5_jobfile.txt --mem-per-cpu 70G --mail-type FAIL,END
# dsq --job-file phastBias_chr6_jobfile.txt --mem-per-cpu 70G --mail-type FAIL,END
# dsq --job-file phastBias_chr7_jobfile.txt --mem-per-cpu 70G --mail-type FAIL,END
# dsq --job-file phastBias_chr8_jobfile.txt --mem-per-cpu 60G --mail-type FAIL,END
# dsq --job-file phastBias_chr9_jobfile.txt --mem-per-cpu 50G --mail-type FAIL,END
# dsq --job-file phastBias_chr10_jobfile.txt --mem-per-cpu 50G --mail-type FAIL,END
# dsq --job-file phastBias_chr11_jobfile.txt --mem-per-cpu 50G --mail-type FAIL,END
# dsq --job-file phastBias_chr12_jobfile.txt --mem-per-cpu 50G --mail-type FAIL,END
# dsq --job-file phastBias_chr13_jobfile.txt --mem-per-cpu 40G --mail-type FAIL,END
# dsq --job-file phastBias_chr14_jobfile.txt --mem-per-cpu 40G --mail-type FAIL,END
# dsq --job-file phastBias_chr15_jobfile.txt --mem-per-cpu 40G --mail-type FAIL,END
# dsq --job-file phastBias_chr16_jobfile.txt --mem-per-cpu 40G --mail-type FAIL,END
# dsq --job-file phastBias_chr17_jobfile.txt --mem-per-cpu 40G --mail-type FAIL,END
# dsq --job-file phastBias_chr18_jobfile.txt --mem-per-cpu 30G --mail-type FAIL,END
# dsq --job-file phastBias_chr19_jobfile.txt --mem-per-cpu 20G --mail-type FAIL,END
# dsq --job-file phastBias_chr20_jobfile.txt --mem-per-cpu 30G --mail-type FAIL,END
# dsq --job-file phastBias_chr21_jobfile.txt --mem-per-cpu 20G --mail-type FAIL,END
# dsq --job-file phastBias_chr22_jobfile.txt --mem-per-cpu 20G --mail-type FAIL,END
# dsq --job-file phastBias_chrX_jobfile.txt --mem-per-cpu 50G --mail-type FAIL,END
