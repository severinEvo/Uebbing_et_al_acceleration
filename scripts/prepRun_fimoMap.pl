#!/usr/bin/perl
use strict;
use warnings;
use utf8;

my @genomes = ("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
	"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
	"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1",
	"tarSyr2","otoGar3","micMur3","proCoq1","galVar1","tupChi1","jacJac1","micOch1",
	"criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6",
	"HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3",
	"chiLan1","octDeg1","speTri2","HLmarMar1","oryCun2","ochPri3","vicPac2",
	"HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1","HLdelLeu1",
	"lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1",
	"bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1",
	"HLcerEla1","susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8",
	"HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1",
	"HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1",
	"HLmanJav1","pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1",
	"myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1","eriEur2","sorAra2",
	"conCri1","loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1",
	"oryAfe1","dasNov3","HLchoHof2","monDom5","sarHar1","HLphaCin1","ornAna2");

foreach my $g (@genomes){
	my @jolma;
	my @a = (1..6);
	foreach my $i (@a){
		my $infile = 'Jolma_' . $i . '.tsv';
		open(IN, "<$infile");
		while(<IN>){
			chomp;
			my @t = split /\t/;
			push @jolma, $t[3];
		}
		close(IN);
		open(OUT, ">>fimo_job_$g\.txt");
		foreach my $j (@jolma){
			print OUT ". ~/.bashrc; mkdir -p data/meme/$g ; cd data/meme/$g ; ";
			print OUT "module load MEME; ";
			print OUT "fimo --bfile ../genomes/$g\_markov-bg-model.txt --motif $j --max-stored-scores 10000000 --oc $j\_out annotations/meme-database/table_s3-$i\.meme $g\.fa; ";
			print OUT "rm $j\_out/*.xml $j\_out/*.tsv $j\_out/*.html\n";
		}
		close(OUT);
		@jolma = ();
	}
	
	my @Hhuman;
	open(IN, "<HOCOMOCO_HUMAN.tsv");
	while(<IN>){
		chomp;
		my @t = split /\t/;
		push @Hhuman, $t[3];
	}
	close(IN);
	my @Hmouse;
	open(IN, "<HOCOMOCO_MOUSE.tsv");
	while(<IN>){
		chomp;
		my @t = split /\t/;
		push @Hmouse, $t[3];
	}
	close(IN);
	my @jaspar;
	open(IN, "<jaspar_modified.tsv");
	while(<IN>){
		chomp;
		my @t = split /\t/;
		my @u = split /_/, $t[3];
		push @jaspar, $u[1];
	}
	close(IN);
	
	open(OUT, ">>fimo_job_$g\.txt");
	foreach my $i (@Hmouse){
		print OUT ". ~/.bashrc; mkdir -p data/meme/$g ; cd data/meme/$g ; ";
		print OUT "module load MEME; ";
		print OUT "fimo --bfile ../genomes/$g\_markov-bg-model.txt --motif $i --max-stored-scores 10000000 --oc $i\_out annotations/meme-database/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme $g\.fa; ";
		print OUT "rm $i\_out/*.xml $i\_out/*.tsv $i\_out/*.html\n";
	}
	close(OUT);	
	
	open(OUT, ">>fimo_job_$g\.txt");
	foreach my $i (@Hhuman){
		print OUT ". ~/.bashrc; mkdir -p data/meme/$g ; cd data/meme/$g ; ";
		print OUT "module load MEME; ";
		print OUT "fimo --bfile ../genomes/$g\_markov-bg-model.txt --motif $i --max-stored-scores 10000000 --oc $i\_out annotations/meme-database/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme $g\.fa; ";
		print OUT "rm $i\_out/*.xml $i\_out/*.tsv $i\_out/*.html\n";
	}
	close(OUT);
	
	open(OUT, ">>fimo_job_$g\.txt");
	foreach my $i (@jaspar){
		print OUT ". ~/.bashrc; mkdir -p data/meme/$g ; cd data/meme/$g ; ";
		print OUT "module load MEME; ";
		print OUT "fimo --bfile ../genomes/$g\_markov-bg-model.txt --motif $i --max-stored-scores 10000000 --oc $i\_out annotations/meme-database/$i\.meme $g\.fa; ";
		print OUT "rm $i\_out/*.xml $i\_out/*.tsv $i\_out/*.html\n";
	}
	close(OUT);
}

# Create jobs:
# module load dSQ
# for i in $(ls fimo_job_*.txt); do echo $i; dsq --job-file $i --mem-per-cpu 300M --mail-type FAIL,END; done
