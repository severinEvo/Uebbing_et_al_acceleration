setwd("data/")

reqPackages<-c("biomaRt")
newPackages<-reqPackages[!(reqPackages %in% utils::installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(newPackages)

# load data
datq<-read.delim("phyloP_out/phyloP-all.tsv.gz")
datq$pos<-paste(datq$chr,":",datq$start,"-",datq$end,sep="")
sign<-datq[p.adjust(datq$pval,method="BH")<.05,]

# prepare table for Fig. 2A
allTab<-data.frame(table(datq$pos))
colnames(allTab)<-c("pos","N")
allTab$pos<-as.character(allTab$pos)
signTab<-data.frame(table(sign$pos))
colnames(signTab)<-c("pos","sign")
signTab$pos<-as.character(signTab$pos)
relAcc<-merge(allTab,signTab,all.x=T)
relAcc$sign<-replace(relAcc$sign,is.na(relAcc$sign),0)
relAcc$relAcc<-relAcc$sign/relAcc$N
pC<-unique(datq[,c("phastCons","pos")])
relAcc<-merge(relAcc,pC)

# correlation test for Fig. 2A
cor.test(relAcc$phastCons,relAcc$relAcc,method="spearman")
# rho = -0.1364883, S = 1.551e+16, p-value = 0

# bin data into 2% bins and plot
list<-split(relAcc$phastCons,cut(relAcc$relAcc,seq(0,.1,.02)))
pdf("../relAcc-vs-phastCons-boxplot.pdf") # Fig. 2A
plot(1,type="n",ylim=0:1,xlim=0:1,las=1,ylab="phastCons",xlab="Acceleration per PRE (%)",
	xaxt="n",bty="n")
axis(1,at=seq(.1,.9,.2),labels=c("0-2","2-4","4-6","6-8","8-10"),tick=F)
segments(y0=unlist(lapply(list,quantile,.05)),y1=unlist(lapply(list,quantile,.95)),
	x0=seq(.1,.9,.2),lwd=2)
segments(y0=unlist(lapply(list,median)),x0=seq(.06,.86,.2),x1=seq(.14,.94,.2),lwd=3)
rect(ybottom=unlist(lapply(list,quantile,.25)),ytop=unlist(lapply(list,quantile,.75)),
	xleft=seq(.06,.86,.2),xright=seq(.14,.94,.2),lwd=2)
legend("topright",c(expression(paste(rho," = -0.14")),expression("P < 10"^-307)),bty="n")
dev.off()

# load connectivity annotations
HiC<-read.delim("PRE-gene-dataset/connection-collection.bed.gz",header=F)
colnames(HiC)<-c("chr","start","end","ensembl_gene_id","PRE_type","connection")
HiC$pos<-paste(HiC$chr,":",HiC$start,"-",HiC$end,sep="")
HiC<-merge(HiC,sign,all.x=T)

# gene table with PRE acceleration stats
all_ubranch<-unique(HiC[,c("pos","branch")])
dummy<-all_ubranch[is.na(all_ubranch$branch),]
all_ubranch<-unique(HiC[!is.na(HiC$branch),c("pos","branch")])
all_repPRE<-data.frame(table(all_ubranch$pos))
colnames(all_repPRE)<-colnames(dummy)<-c("pos","nAccBranches")
all_repPRE$pos<-as.character(all_repPRE$pos)
dummy$nAccBranches<-0
all_repPRE<-rbind(all_repPRE,dummy)
all_ugene<-unique(HiC[,c("pos","ensembl_gene_id")])
all_nPREs<-data.frame(table(all_ugene$ensembl_gene_id))
colnames(all_nPREs)<-c("ensembl_gene_id","nPREs")
all_nPREs$ensembl_gene_id<-as.character(all_nPREs$ensembl_gene_id)
all<-merge(unique(HiC[,c("pos","ensembl_gene_id")]),all_nPREs)
all<-merge(all,all_repPRE)
dummy<-data.frame(table(all$ensembl_gene_id[all$nAccBranches!=0]))
colnames(dummy)<-c("ensembl_gene_id","nAccPREs")
all<-merge(all,dummy,all.x=T)
all$nAccPREs<-replace(all$nAccPREs,is.na(all$nAccPREs),0)
genePREtab<-unique(all[,c("ensembl_gene_id","nPREs","nAccPREs")])
genePREtab$relAcc<-genePREtab$nAccPREs/genePREtab$nPREs

# correlation test on raw data for Fig. 2B
cor.test(genePREtab$nPREs,genePREtab$relAcc,method="spearman")
# S = 7.2855e+12, p-value = 0, rho = 0.2810984

# bin data for plotting Fig. 2B
nTable<-data.frame(nPREs=1:200,percAcc=double(200))
dummy<-data.frame(table(genePREtab$nPREs))
colnames(dummy)<-c("nPREs","n")
dummy$nPREs<-as.integer(dummy$nPREs)
nTable<-merge(nTable,dummy,all.x=T)
for(i in 1:200){
	nTable$percAcc[i]<-sum(genePREtab$nAccPREs[genePREtab$nPREs==i])/
	sum(genePREtab$nPREs[genePREtab$nPREs==i])
}

# regression line for Fig. 2B
sdy<-sd(genePREtab$relAcc,na.rm=T)
sdx<-sd(genePREtab$nPREs)
cc<-cor.test(genePREtab$nPREs,genePREtab$relAcc)$est
b1<-cc*(sdy/sdx)
meany<-mean(genePREtab$relAcc,na.rm=T)
meanx<-mean(genePREtab$nPREs)
b0<-meany-b1*meanx

pdf("../nPREs-v-relAcc.pdf") # Fig. 2B
plot(nTable$nPREs,nTable$percAcc,xlim=c(0,99),ylim=c(.15,.35),las=1,bty="n",
	xlab="PREs per gene",ylab="Accelerated PREs (%)",yaxt="n",pch=16)
axis(2,at=seq(.15,.35,.05),labels=seq(15,35,5),las=1)
abline(b0,b1,lwd=2)
legend("topleft",c(expression(paste(rho," = 0.281")),expression("P < 10"^-307)),bty="n")
dev.off()

if(file.exists("PRE-gene-dataset/geneSignTable.tsv")){
	geneSign<-read.delim("PRE-gene-dataset/geneSignTable.tsv")
}else{
	uHiC<-unique(HiC[,c("pos","ensembl_gene_id","branch")])
	uHiCtab<-merge(uHiC,datq[,c("branch","pos")])
	uHiCsign<-merge(uHiC,sign[c("branch","pos")])
	geneTest<-data.frame(table(uHiCtab$ensembl_gene_id))
	colnames(geneTest)<-c("ensembl_gene_id","tested")
	geneTest$ensembl_gene_id<-as.character(geneTest$ensembl_gene_id)
	geneSign<-data.frame(table(uHiCsign$ensembl_gene_id))
	colnames(geneSign)<-c("ensembl_gene_id","sign")
	geneSign$ensembl_gene_id<-as.character(geneSign$ensembl_gene_id)
	geneSign<-merge(geneTest,geneSign,all.x=T)
	geneSign$sign<-replace(geneSign$sign,is.na(geneSign$sign),0)
	geneSign$relAcc<-geneSign$sign/geneSign$tested
	write.table(geneSign,"PRE-gene-dataset/geneSignTable.tsv",row.names=F,sep="\t",quote=F)	
}
all<-merge(all,geneSign[,-(2:3)])

# load and integrate gnomAD data
gnomad<-read.delim("../annotations/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz")
colnames(gnomad)[2]<-"ensembl_transcript_id"
library(biomaRt)
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
transl<-getBM(filters=c("ensembl_transcript_id"),values=gnomad$ensembl_transcript_id,
	attributes=c("ensembl_gene_id","ensembl_transcript_id"),mart=ensembl)
gnomad<-merge(gnomad[,c("gene","ensembl_transcript_id","oe_mis","oe_mis_pphen","oe_lof")],transl)
gnomAData<-merge(all,gnomad)

cor.test(gnomAData$nAccBranches,gnomAData$oe_mis_pphen,method="spearman")
cor.test(gnomAData$nAccBranches,gnomAData$oe_lof,method="spearman")
# slight positive correlations: more constrained genes have fewer acceleration events

# correlation tests for Fig. 2C
cor.test(gnomAData$nPREs,gnomAData$oe_mis_pphen,method="spearman")
cor.test(gnomAData$nPREs,gnomAData$oe_lof,method="spearman")
# negative correlations: more constrained genes have more enhancers

cor.test(gnomAData$nAccPREs,gnomAData$oe_mis_pphen,method="spearman")
cor.test(gnomAData$nAccPREs,gnomAData$oe_lof,method="spearman")
# negative correlations: more constrained genes have more accelerated enhancers
# this is a combination effect of the two effects above

# correlations tests for Fig. 2D
cor.test(gnomAData$relAcc,gnomAData$oe_mis_pphen,method="spearman")
cor.test(gnomAData$relAcc,gnomAData$oe_lof,method="spearman")
# slight positive correlations: more constrained genes show less acceleration

pdf("../../gnomAD.pdf") # Fig. 2C,D
par(mfrow=c(2,2),las=1,mar=c(5,3,3,0))
list1<-split(gnomAData$oe_mis_pphen[!is.na(gnomAData$oe_mis_pphen)],
	cut(gnomAData$nPREs[!is.na(gnomAData$oe_mis_pphen)],seq(0,150,5)))
plot(1,1,type="n",xlim=c(1,30),ylim=c(0,1.5),xaxt="n",yaxt="n",bty="n",
	ylab="",xlab="CREs per gene",main="PolyPhen-2 missense variants")
axis(1,c(1,50,100,150),at=c(1,10,20,30))
polygon(x=c(1:30,30:1),border=NA,col=gray(.8),
	y=c(unlist(lapply(list1,quantile,.025)),unlist(lapply(list1,quantile,.975))[30:1]))
polygon(x=c(1:30,30:1),border=NA,col=gray(.6),
	y=c(unlist(lapply(list1,quantile,.25)),unlist(lapply(list1,quantile,.75))[30:1]))
lines(unlist(lapply(list1,median)),type="l",lwd=2)
legend("topleft",c("rho = -0.17","P = 0"),bty="n")
list2<-split(gnomAData$oe_lof[!is.na(gnomAData$oe_lof)],
	cut(gnomAData$nPREs[!is.na(gnomAData$oe_lof)],seq(0,150,5)))
plot(1,1,type="n",xlim=c(1,30),ylim=c(0,1.5),xaxt="n",yaxt="n",bty="n",
	ylab="",xlab="CREs per gene",main="Predicted loss-of-function variants")
axis(1,c(1,50,100,150),at=c(1,10,20,30))
polygon(x=c(1:30,30:1),border=NA,col=gray(.8),
	y=c(unlist(lapply(list2,quantile,.025)),unlist(lapply(list2,quantile,.975))[30:1]))
polygon(x=c(1:30,30:1),border=NA,col=gray(.6),
	y=c(unlist(lapply(list2,quantile,.25)),unlist(lapply(list2,quantile,.75))[30:1]))
lines(unlist(lapply(list2,median)),type="l",lwd=2)
legend("topleft",c("rho = -0.26","P = 0"),bty="n")
list1<-split(gnomAData$oe_mis_pphen[!is.na(gnomAData$oe_mis_pphen)],
	cut(gnomAData$relAcc[!is.na(gnomAData$oe_mis_pphen)],(0:10)/200))
plot(1,1,type="n",xlim=c(1,10),ylim=c(0,1.5),xaxt="n",yaxt="n",bty="n",
	ylab="",xlab="Relative acceleration per gene")
axis(1,at=1:10,labels=(1:10)/200)
polygon(x=c(1:10,10:1),border=NA,col=gray(.8),
	y=c(unlist(lapply(list1,quantile,.025)),unlist(lapply(list1,quantile,.975))[10:1]))
polygon(x=c(1:10,10:1),border=NA,col=gray(.6),
	y=c(unlist(lapply(list1,quantile,.25)),unlist(lapply(list1,quantile,.75))[10:1]))
lines(unlist(lapply(list1,median)),type="l",lwd=2)
legend("topleft",c("rho = 0.054","P = 3.4e-180"),bty="n")
list2<-split(gnomAData$oe_lof[!is.na(gnomAData$oe_lof)],
	cut(gnomAData$relAcc[!is.na(gnomAData$oe_lof)],(0:10)/200))
plot(1,1,type="n",xlim=c(1,10),ylim=c(0,1.5),xaxt="n",yaxt="n",bty="n",
	ylab="",xlab="Relative acceleration per gene")
axis(1,at=1:10,labels=(1:10)/200)
polygon(x=c(1:10,10:1),border=NA,col=gray(.8),
	y=c(unlist(lapply(list2,quantile,.025)),unlist(lapply(list2,quantile,.975))[10:1]))
polygon(x=c(1:10,10:1),border=NA,col=gray(.6),
	y=c(unlist(lapply(list2,quantile,.25)),unlist(lapply(list2,quantile,.75))[10:1]))
lines(unlist(lapply(list2,median)),type="l",lwd=2)
legend("topleft",c("rho = 0.051","P = 6.4e-159"),bty="n")
dev.off()

q(save="no")
