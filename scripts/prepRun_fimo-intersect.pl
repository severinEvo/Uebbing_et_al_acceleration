#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @chr = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX");
my @spec = ("panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5",
	"macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1","colAng1",
	"HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2","otoGar3","micMur3",
	"proCoq1","galVar1","tupChi1","jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10",
	"HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2",
	"HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2","HLmarMar1","oryCun2","ochPri3",
	"vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1","HLdelLeu1",
	"lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1","bosMut1",
	"bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1","HLcerEla1",
	"susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1","panTig1",
	"HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1",
	"odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1","pteAle1","HLpteVam2","rouAeg1",
	"HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1",
	"eriEur2","sorAra2","conCri1","loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2",
	"eleEdw1","oryAfe1","dasNov3","HLchoHof2","monDom5","sarHar1","HLphaCin1","ornAna2");

foreach my $c (@chr){
	open(OUT, ">$c\_fimo-intersect_job.txt");
	print OUT ". ~/.bashrc; cd data/meme/$c; ";
	print OUT "module load BEDTools; ";
	print OUT "bedtools intersect -f 1 -b ../../MAFs/tmp_mafs/$c\/PREs_$c\.bed -a ../hg38/fimo_hg38-all.bed.gz >fimo_hg38-$c\_intersect.bed; ";
	print OUT "LIST=`cat phyloP-$c\_all.lst`; for i in \$LIST; do ";
	print OUT "bedtools intersect -f 1 -wo -a fimo_hg38-$c\_intersect.bed -b ../../MAFs/tmp_mafs/$c\/\$i\\.bed | cut -f 1-6,10 >>\$i\\_fimo_all_$c\.bed; ";
	print OUT "done\n";
	foreach my $s (@spec){
		print OUT ". ~/.bashrc; cd data/meme/$c; ";
		print OUT "grep $s $c\:*_trans.bed | cut -f 3 -d ':' | sort -k1,1 -k2,2n >$s\_$c\_trans.bed; ";
		print OUT "module load BEDTools; ";
		print OUT "bedtools intersect -f 1 -b $s\_$c\_trans.bed -a ../$s\/fimo_$s\-all.bed.gz >fimo_$s\-$c\_intersect.bed; ";
		print OUT "LIST=`cat phyloP-$c\_all.lst`; for i in \$LIST; do ";
		print OUT "grep $s \$i\\_trans.bed | bedtools intersect -f 1 -wo -a fimo_$s\-$c\_intersect.bed -b - | cut -f 1-6,10 >>\$i\\_fimo_all_$c\.bed; ";
		print OUT "done\n";
	}
	close(OUT);
}

# for i in $(ls chr*_fimo-intersect_job.txt); do dsq --jobfile $i --mem-per-cpu 100M --mail-type FAIL,END; done
# for i in $(ls dsq-chr*.sh); do sbatch $i ; done
