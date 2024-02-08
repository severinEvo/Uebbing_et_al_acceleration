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

open(OUT, ">download-prep-genomes.txt");
foreach my $g (@genomes){
	print OUT ". ~/.bashrc; mkdir -p data/genomes; cd data/genomes/; ";
	print OUT "wget https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/Genome2bitFiles/$g\.2bit; ";
	print OUT "twoBitToFa $g\.2bit $g\.fa; ";
	print OUT "module load MEME; ";
	print OUT "fasta-get-markov $g\.fa $g\_markov-bg-model.txt\n";
close(OUT);

# Create jobs:
# module load dSQ
# dsq --jobfile download-prep-genomes.txt --mem-per-cpu 100M --mail-type FAIL,END
