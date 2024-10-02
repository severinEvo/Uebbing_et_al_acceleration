setwd("data/")

reqPackages<-c("BiocManager","gam")
reqBcPackages<-c("biomaRt")
newPackages<-reqPackages[!(reqPackages %in% utils::installed.packages()[,"Package"])]
if(length(newPackages)>0) install.packages(newPackages)
newBcPackages<-reqBcPackages[!(reqBcPackages %in% utils::installed.packages()[,"Package"])]
if(length(newBcPackages)>0) BiocManager::install(newBcPackages)

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
dummy$ensembl_gene_id<-as.character(dummy$ensembl_gene_id)
all<-merge(all,dummy,all.x=T)
all$nAccPREs<-replace(all$nAccPREs,is.na(all$nAccPREs),0)

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

tau<-read.delim("GTEx_tau.tsv")
gnomAData<-merge(gnomAData,tau)

# tests for Fig. 2B-D
library(gam)
gamRelAcc<-gam(relAcc ~ lo(nPREs)+lo(oe_mis_pphen)+lo(oe_lof)+lo(tau),data=gnomAData)
summary(gamRelAcc)
statsRelAcc<-summary(gamRelAcc)
statsRelAcc$anova$'Pr(F)'

gamN<-gam(nPREs ~ lo(oe_mis_pphen)+lo(oe_lof)+lo(tau),data=gnomAData)
summary(gamN)
statsN<-summary(gamN)
statsN$anova$'Pr(F)'

# plot Fig. 2B-D
pdf("../GAM.pdf",width=12,height=16)
par(mfrow=c(4,3),las=1)
plot(1,type="n")
plot(density(gnomAData$nPREs))
plot(gamRelAcc,se=T,rugplot=F)
plot(gamN,se=T,rugplot=F)
plot(density(gnomAData$oe_mis_pphen,na.rm=T))
plot(density(gnomAData$oe_lof,na.rm=T))
plot(density(gnomAData$tau,na.rm=T))
dev.off()

q(save="no")
