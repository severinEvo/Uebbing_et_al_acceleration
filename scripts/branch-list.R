pan<-list("panTro5","panPan2","hg38")
hominini<-list("hg38",c("panTro5","panPan2"),"gorGor5")
homininae<-list(c("hg38","panTro5","panPan2"),"gorGor5","ponAbe2")
hominidae<-list(c("hg38","panTro5","panPan2","gorGor5"),"ponAbe2","nomLeu3")
hominoidea<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2"),"nomLeu3",
	c("rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1",
		"rhiRox1","rhiBie1","colAng1","HLpilTep1"))
papionini<-list(c("rheMac8","macFas5","macNem1"),c("papAnu3","manLeu1","cerAty1"),"chlSab2")
cercopithecinae<-list(c("rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1"),"chlSab2",
	c("nasLar1","rhiRox1","rhiBie1","colAng1","HLpilTep1"))
presbytini<-list("nasLar1",c("rhiRox1","rhiBie1"),c("colAng1","HLpilTep1"))
colobini<-list("colAng1","HLpilTep1",c("nasLar1","rhiRox1","rhiBie1"))
colobinae<-list(c("nasLar1","rhiRox1","rhiBie1"),c("colAng1","HLpilTep1"),
	c("rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2"))
cercopithecidae<-list(c("rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2"),
	c("nasLar1","rhiRox1","rhiBie1","colAng1","HLpilTep1"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3"))
catarrhini<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3"),
	c("rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1",
		"rhiRox1","rhiBie1","colAng1","HLpilTep1"),
	c("calJac3","aotNan1","saiBol1","cebCap1"))
callitrichoidea<-list("calJac3","aotNan1",c("saiBol1","cebCap1"))
cebidae<-list("saiBol1","cebCap1",c("calJac3","aotNan1"))
platyrrhini<-list(c("calJac3","aotNan1"),c("saiBol1","cebCap1"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5","macNem1",
		"papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1","colAng1",
		"HLpilTep1"))
simiiformes<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1"),
	c("calJac3","aotNan1","saiBol1","cebCap1"),"tarSyr2")
haplorrhini<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1"),
	"tarSyr2",c("otoGar3","micMur3","proCoq1"))
lemuroidea<-list("micMur3","proCoq1","otoGar3")
strepsirrhini<-list("otoGar3",c("micMur3","proCoq1"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5","macNem1",
		"papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1","colAng1",
		"HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2"))
primates<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5",
		"macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1",
		"colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2"),
	c("otoGar3","micMur3","proCoq1"),"galVar1")
primatomorpha<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2",
		"otoGar3","micMur3","proCoq1"),
	"galVar1","tupChi1")
euarchonta<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2",
		"otoGar3","micMur3","proCoq1","galVar1"),"tupChi1",
	c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1",
		"rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3",
		"chiLan1","octDeg1","speTri2","HLmarMar1","oryCun2","ochPri3"))
'arvicolinae-cricetinae'<-list("micOch1",c("criGri1","mesAur1"),"perManBai1")
cricetidae<-list(c("micOch1","criGri1","mesAur1"),"perManBai1",
	c("mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1"))
mus<-list(c("mm10","HLmusCar1"),"HLmusPah1","rn6")
murinae<-list(c("mm10","HLmusCar1","HLmusPah1"),"rn6","HLmerUng1")
muridae<-list(c("mm10","HLmusCar1","HLmusPah1","rn6"),"HLmerUng1",
	c("micOch1","criGri1","mesAur1","perManBai1"))
eumuroida<-list(c("micOch1","criGri1","mesAur1","perManBai1"),c("mm10","HLmusCar1",
		"HLmusPah1","rn6","HLmerUng1"),
	"nanGal1")
muroidea<-list(c("micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1",
		"rn6","HLmerUng1"),
	"nanGal1","jacJac1")
myodonta<-list("jacJac1",
	c("micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6",
		"HLmerUng1","nanGal1"),
	c("HLcasCan1","HLdipOrd2"))
castorimorpha<-list("HLcasCan1","HLdipOrd2",
	c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6",
		"HLmerUng1","nanGal1"))
phiomorpha<-list("hetGla2","HLfukDam1",c("cavPor3","chiLan1","octDeg1"))
'octodontoidea-chinchilloidea'<-list("chiLan1","octDeg1","cavPor3")
caviomorpha<-list("cavPor3",c("chiLan1","octDeg1"),c("hetGla2","HLfukDam1"))
hystricomorpha<-list(c("hetGla2","HLfukDam1"),c("cavPor3","chiLan1","octDeg1"),
	c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6",
		"HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2"))
'myodonta-castorimorpha-hystricomorpha'<-list(c("jacJac1","micOch1","criGri1","mesAur1",
	"perManBai1","mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1",
		"HLdipOrd2"),
	c("hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1"),c("speTri2","HLmarMar1"))
sciuromorpha<-list("speTri2","HLmarMar1",
	c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6",
		"HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1",
		"octDeg1"))
rodentia<-list(c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1",
		"HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1",
		"cavPor3","chiLan1","octDeg1"),
	c("speTri2","HLmarMar1"),c("oryCun2","ochPri3"))
lagomorpha<-list("oryCun2","ochPri3",
	c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6",
		"HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1",
		"octDeg1","speTri2","HLmarMar1"))
glires<-list(c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1",
		"HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1",
		"cavPor3","chiLan1","octDeg1","speTri2","HLmarMar1"),
	c("oryCun2","ochPri3"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2",
		"otoGar3","micMur3","proCoq1","galVar1","tupChi1"))
euarchontoglires<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2",
		"otoGar3","micMur3","proCoq1","galVar1","tupChi1"),
	c("jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1",
		"rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3",
		"chiLan1","octDeg1","speTri2","HLmarMar1","oryCun2","ochPri3"),
	c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1","HLdelLeu1",
		"lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1","bosMut1",
		"bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1","HLcerEla1",
		"susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1",
		"panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1",
		"ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1","pteAle1",
		"HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2",
		"HLminNat1","HLdesRot1","eriEur2","sorAra2","conCri1"))
camelus<-list(c("HLcamFer2","HLcamBac1"),"HLcamDro1","vicPac2")
camelidae<-list("vicPac2",c("HLcamFer2","HLcamBac1","HLcamDro1"),
	c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8",
		"HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1",
		"HLodoVir1","HLcerEla1","susScr11"))
delphinidae<-list("HLturTru3","HLorcOrc1","HLdelLeu1")
delphinoidea<-list(c("HLturTru3","HLorcOrc1"),"HLdelLeu1","lipVex1")
'delphinoidea-inoidea'<-list(c("HLturTru3","HLorcOrc1","HLdelLeu1"),"lipVex1","phyCat1")
odontoceti<-list(c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1"),"phyCat1",
	c("balAcu1","HLbalMys1"))
mysticeti<-list("balAcu1","HLbalMys1",c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1"))
cetacea<-list(c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1"),c("balAcu1","HLbalMys1"),
	c("bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2",
		"panHod1","HLodoVir1","HLcerEla1"))
bovina<-list(c("bosTau8","HLbosInd1"),c("bisBis1","bosMut1"),"bubBub1")
bovini<-list(c("bosTau8","HLbosInd1","bisBis1","bosMut1"),"bubBub1",
	c("HLoviAri4","HLoviCan1","HLcapHir2","panHod1"))
ovis<-list("HLoviAri4","HLoviCan1","HLcapHir2")
caprinae<-list(c("HLoviAri4","HLoviCan1"),"HLcapHir2","panHod1")
caprini<-list(c("HLoviAri4","HLoviCan1","HLcapHir2"),"panHod1",
	c("bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1"))
bovidae<-list(c("bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1"),
	c("HLoviAri4","HLoviCan1","HLcapHir2","panHod1"),c("HLodoVir1","HLcerEla1"))
cervidae<-list("HLodoVir1","HLcerEla1",
	c("bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2",
		"panHod1"))
pecora<-list(c("bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1",
		"HLcapHir2","panHod1"),
	c("HLodoVir1","HLcerEla1"),c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1",
		"balAcu1","HLbalMys1"))
cetruminantia<-list(c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1",
		"HLbalMys1"),
	c("bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2",
		"panHod1","HLodoVir1","HLcerEla1"),
	"susScr11")
artiofabula<-list(c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1",
		"HLbalMys1","bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4",
		"HLoviCan1","HLcapHir2","panHod1","HLodoVir1","HLcerEla1"),
	"susScr11",c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1"))
cetartiodactyla<-list(c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1"),
	c("HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8",
		"HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1",
		"HLodoVir1","HLcerEla1","susScr11"),
	c("HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1","panTig1","HLpanPar1",
		"canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1",
		"lepWed1","neoSch1","manPen1","HLmanJav1","pteAle1","HLpteVam2","rouAeg1","HLrhiSin1",
		"HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1"))
equus<-list("HLequCab3","equPrz1","HLequAsi1")
equidae<-list(c("HLequCab3","equPrz1"),"HLequAsi1","cerSim1")
perissodactyla<-list(c("HLequCab3","equPrz1","HLequAsi1"),"cerSim1",
	c("felCat8","HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1",
		"HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1"))
felinae<-list("felCat8","HLaciJub1",c("panTig1","HLpanPar1"))
panthera<-list("panTig1","HLpanPar1",c("felCat8","HLaciJub1"))
felidae<-list(c("felCat8","HLaciJub1"),c("panTig1","HLpanPar1"),
	c("canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1",
		"lepWed1","neoSch1"))
canidae<-list("canFam3","HLlycPic1",
	c("musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1"))
mustelidae<-list("musFur1","HLenhLut1","HLailFul1")
musteloidea<-list(c("musFur1","HLenhLut1"),"HLailFul1",
	c("ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1"))
ursidae<-list("ailMel1","ursMar1",c("odoRosDiv1","lepWed1","neoSch1"))
phocidae<-list("lepWed1","neoSch1","odoRosDiv1")
pinnipedia<-list("odoRosDiv1",c("lepWed1","neoSch1"),c("ailMel1","ursMar1"))
arctoidea<-list(c("musFur1","HLenhLut1","HLailFul1"),c("ailMel1","ursMar1","odoRosDiv1",
		"lepWed1","neoSch1"),
	c("canFam3","HLlycPic1"))
caniformia<-list(c("canFam3","HLlycPic1"),
	c("musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1"),
	c("felCat8","HLaciJub1","panTig1","HLpanPar1"))
carnivora<-list(c("felCat8","HLaciJub1","panTig1","HLpanPar1"),
	c("canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1",
		"lepWed1","neoSch1"),c("manPen1","HLmanJav1"))
pholidota<-list("manPen1","HLmanJav1",
	c("felCat8","HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1",
		"HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1"))
ferae<-list(c("felCat8","HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1",
		"HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1"),
	c("manPen1","HLmanJav1"),c("HLequCab3","equPrz1","HLequAsi1","cerSim1"))
pteropus<-list("pteAle1","HLpteVam2","rouAeg1")
pteropodidae<-list(c("pteAle1","HLpteVam2"),"rouAeg1",c("HLrhiSin1","HLhipArm1"))
yinochiroptera<-list("HLrhiSin1","HLhipArm1",c("pteAle1","HLpteVam2","rouAeg1"))
yinpterochiroptera<-list(c("pteAle1","HLpteVam2","rouAeg1"),c("HLrhiSin1","HLhipArm1"),
	c("eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1"))
vespertilionoidea<-list(c("eptFus1","myoDav1","myoBra1","myoLuc2"),"HLminNat1","HLdesRot1")
yangochiroptera<-list(c("eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1"),"HLdesRot1",
	c("pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1"))
chiroptera<-list(c("pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1"),
	c("eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1"),
	c("HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1","panTig1","HLpanPar1",
		"canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1",
		"odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1"))
scrotifera<-list(c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1",
		"HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1",
		"bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1",
		"HLcerEla1","susScr11"),
	c("HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1","panTig1","HLpanPar1",
		"canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1",
		"odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1","pteAle1","HLpteVam2",
		"rouAeg1","HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1",
		"HLdesRot1"),c("eriEur2","sorAra2","conCri1"))
'erinaceidae-soricidae'<-list("eriEur2","sorAra2","conCri1")
eulipotyphla<-list(c("eriEur2","sorAra2"),"conCri1",
	c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1","HLdelLeu1",
		"lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1","bosMut1",
		"bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1","HLcerEla1",
		"susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1",
		"panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1",
		"ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1","pteAle1",
		"HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2",
		"HLminNat1","HLdesRot1"))
laurasiatheria<-list(c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1",
		"HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1",
		"bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1",
		"HLcerEla1","susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8",
		"HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1",
		"HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1",
		"HLmanJav1","pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1",
		"myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1"),
	c("eriEur2","sorAra2","conCri1"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2",
		"otoGar3","micMur3","proCoq1","galVar1","tupChi1","jacJac1","micOch1","criGri1",
		"mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1",
		"HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2",
		"HLmarMar1","oryCun2","ochPri3"))
boreoeutheria<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8",
		"macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1",
		"rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2",
		"otoGar3","micMur3","proCoq1","galVar1","tupChi1","jacJac1","micOch1","criGri1",
		"mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1",
		"HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2",
		"HLmarMar1","oryCun2","ochPri3"),
	c("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1","HLdelLeu1",
		"lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1","bosMut1",
		"bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1","HLcerEla1",
		"susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8","HLaciJub1",
		"panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1",
		"ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1","pteAle1",
		"HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2",
		"HLminNat1","HLdesRot1","eriEur2","sorAra2","conCri1"),
	c("loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1","oryAfe1","dasNov3",
		"HLchoHof2"))
tethytheria<-list("loxAfr3","triMan1","HLproCap2")
paenungulata<-list(c("loxAfr3","triMan1"),"HLproCap2",c("chrAsi1","echTel2","eleEdw1","oryAfe1"))
afrosoricida<-list("chrAsi1","echTel2","eleEdw1")
afroinsectivora<-list(c("chrAsi1","echTel2"),"eleEdw1","oryAfe1")
afroinsectiphilia<-list(c("chrAsi1","echTel2","eleEdw1"),"oryAfe1",c("loxAfr3","triMan1",
	"HLproCap2"))
afrotheria<-list(c("loxAfr3","triMan1","HLproCap2"),c("chrAsi1","echTel2","eleEdw1","oryAfe1"),
	c("dasNov3","HLchoHof2"))
xenarthra<-list("dasNov3","HLchoHof2",
	c("loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1","oryAfe1"))
atlantogenata<-list(c("loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1","oryAfe1"),
	c("dasNov3","HLchoHof2"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5","macNem1",
		"papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1","colAng1",
		"HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2","otoGar3","micMur3",
		"proCoq1","galVar1","tupChi1","jacJac1","micOch1","criGri1","mesAur1","perManBai1",
		"mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2",
		"hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2","HLmarMar1","oryCun2",
		"ochPri3","vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1",
		"HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1",
		"bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1",
		"HLcerEla1","susScr11","HLequCab3","equPrz1","HLequAsi1","cerSim1","felCat8",
		"HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1",
		"HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1",
		"HLmanJav1","pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1",
		"myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1","eriEur2","sorAra2","conCri1"))
eutheria<-list(c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5",
		"macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1",
		"colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2","otoGar3",
		"micMur3","proCoq1","galVar1","tupChi1","jacJac1","micOch1","criGri1","mesAur1",
		"perManBai1","mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1",
		"HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2","HLmarMar1",
		"oryCun2","ochPri3","vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3",
		"HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8",
		"HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2",
		"panHod1","HLodoVir1","HLcerEla1","susScr11","HLequCab3","equPrz1","HLequAsi1",
		"cerSim1","felCat8","HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1",
		"HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1",
		"manPen1","HLmanJav1","pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1",
		"eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1","eriEur2","sorAra2",
		"conCri1"),
	c("loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1","oryAfe1","dasNov3",
		"HLchoHof2"),
	c("monDom5","sarHar1","HLphaCin1"))
australidelphia<-list("sarHar1","HLphaCin1","monDom5")
marsupialia<-list("monDom5",c("sarHar1","HLphaCin1"),
	c("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5",
		"macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1",
		"colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2","otoGar3",
		"micMur3","proCoq1","galVar1","tupChi1","jacJac1","micOch1","criGri1","mesAur1",
		"perManBai1","mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1",
		"HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2","HLmarMar1",
		"oryCun2","ochPri3","vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3",
		"HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8",
		"HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2",
		"panHod1","HLodoVir1","HLcerEla1","susScr11","HLequCab3","equPrz1","HLequAsi1",
		"cerSim1","felCat8","HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1",
		"HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1",
		"manPen1","HLmanJav1","pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1",
		"eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1","eriEur2","sorAra2",
		"conCri1","loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1","oryAfe1",
		"dasNov3","HLchoHof2"))

branchlist<-list(pan,hominini,homininae,hominidae,hominoidea,papionini,cercopithecinae,
	presbytini,colobini,colobinae,cercopithecidae,catarrhini,callitrichoidea,cebidae,
	platyrrhini,simiiformes,haplorrhini,lemuroidea,strepsirrhini,primates,primatomorpha,
	euarchonta,get('arvicolinae-cricetinae'),cricetidae,mus,murinae,muridae,eumuroida,
	muroidea,myodonta,castorimorpha,phiomorpha,get("octodontoidea-chinchilloidea"),
	caviomorpha,hystricomorpha,get("myodonta-castorimorpha-hystricomorpha"),sciuromorpha,
	rodentia,lagomorpha,glires,euarchontoglires,camelus,camelidae,delphinidae,delphinoidea,
	get("delphinoidea-inoidea"),odontoceti,mysticeti,cetacea,bovina,bovini,ovis,caprinae,
	caprini,bovidae,cervidae,pecora,cetruminantia,artiofabula,cetartiodactyla,equus,equidae,
	perissodactyla,felinae,panthera,felidae,canidae,mustelidae,musteloidea,ursidae,phocidae,
	pinnipedia,arctoidea,caniformia,carnivora,pholidota,ferae,pteropus,pteropodidae,
	yinochiroptera,yinpterochiroptera,vespertilionoidea,yangochiroptera,chiroptera,
	get("erinaceidae-soricidae"),scrotifera,eulipotyphla,laurasiatheria,boreoeutheria,
	tethytheria,paenungulata,afrosoricida,afroinsectivora,afroinsectiphilia,afrotheria,
	xenarthra,atlantogenata,eutheria,australidelphia,marsupialia)
names(branchlist)<-c("pan","hominini","homininae","hominidae","hominoidea","papionini",
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
	"yinpterochiroptera","vespertilionoidea","yangochiroptera","chiroptera",
	"erinaceidae-soricidae","scrotifera","eulipotyphla","laurasiatheria","boreoeutheria",
	"tethytheria","paenungulata","afrosoricida","afroinsectivora","afroinsectiphilia",
	"afrotheria","xenarthra","atlantogenata","eutheria","australidelphia","marsupialia")

save(branchlist,file="../annotations/branch-list.Rdata")

q(save="no")
