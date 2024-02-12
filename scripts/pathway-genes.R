HiC<-read.delim("data/PRE-gene-database/connection-collection.bed.gz",header=F)
colnames(HiC)<-c("chr","start","end","ensembl_gene_id","CRE_type","connection")
HiC$pos<-paste(HiC$chr,":",HiC$start,"-",HiC$end,sep="")
sign<-read.delim("data/phyloP_out/phyloP-sign_0based.tsv.gz")
sign$pos<-paste(sign$chr,":",sign$start,"-",sign$end,sep="")
uHiC<-unique(HiC[,c(4,5,7)])

if(file.exists("data/allGenes.tsv.gz")){
	allGenes<-read.delim("data/allGenes.tsv.gz")
}else{
	library(biomaRt)
	ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
	allGenes<-getBM(filters=c("ensembl_gene_id"),mart=ensembl,
		attributes=c("ensembl_gene_id","external_gene_name"),values=unique(HiC$ensembl_gene_id))
	write.table(allGenes,gzfile("data/allGenes.tsv.gz"),sep="\t",quote=F,row.names=F)
}

lst<-read.delim("annotations/pathways.tsv")
lst<-merge(lst,allGenes)
lst<-merge(lst,uHiC)

datq<-read.delim("phyloP_out/phyloP-all.tsv.gz")
datq$pos<-paste(datq$chr,":",datq$start,"-",datq$end,sep="")
lst<-merge(lst,datq[,c(4,9)])
sign$sign<-1
lst<-merge(lst,sign[,c(4,6,7)],all.x=T)
lst$sign<-replace(lst$sign,is.na(lst$sign),0)

load("annotations/branch-list.Rdata")
branches<-names(branchlist)
pathways<-sort(unique(lst$pathway))

datsign<-merge(uHiC,datq[,c(4,9)])
datsign<-merge(datsign,sign[,c(4,6,7)],all.x=T)
datsign$sign<-replace(datsign$sign,is.na(datsign$sign),0)
geneSignTab<-data.frame(table(datsign[datsign$sign==1,2:3]))
colnames(geneSignTab)[3]<-"sign"
geneSignTab$branch<-as.character(geneSignTab$branch)
geneSignTab$ensembl_gene_id<-as.character(geneSignTab$ensembl_gene_id)
geneTestTab<-data.frame(table(datsign[,2:3]))
colnames(geneTestTab)[3]<-"tested"
geneTestTab$branch<-as.character(geneTestTab$branch)
geneTestTab$ensembl_gene_id<-as.character(geneTestTab$ensembl_gene_id)
geneTab<-merge(geneTestTab,geneSignTab,all=T)
geneTab$sign<-replace(geneTab$sign,is.na(geneTab$sign),0)
geneTab$relAcc<-geneTab$sign/geneTab$tested
geneTab$relAcc<-replace(geneTab$relAcc,is.na(geneTab$relAcc),0)

n<-10000
genePREtab<-data.frame(table(uHiC$ensembl_gene_id))
colnames(genePREtab)<-c("ensembl_gene_id","PREs")
genePREtab$ensembl_gene_id<-as.character(genePREtab$ensembl_gene_id)
uGeneTab<-data.frame(table(datsign$ensembl_gene_id[ datsign$sign==1]))
colnames(uGeneTab)<-c("ensembl_gene_id","accEvents")
uGeneTab$ensembl_gene_id<-as.character(uGeneTab$ensembl_gene_id)
if(file.exists("data/resMean.Rdata")){
	load("data/resMean.Rdata")
}else{
	genes<-unique(geneTab$ensembl_gene_id)
	resData<-data.frame(meanRelAcc=double(n),meanNoPREs=double(n),meanNoAcc=double(n))
	for(i in 1:n){
		if(i %% 100==0){
			print(i)
		}
		dummy<-sample(genes,50)
		dummy2<-geneTab[geneTab$ensembl_gene_id %in% dummy,]
		resData[i,1]<-mean(dummy2$relAcc,na.rm=T)
		dummy2<-genePREtab[genePREtab$ensembl_gene_id %in% dummy,]
		resData[i,2]<-mean(dummy2$PREs)
		dummy2<-uGeneTab[uGeneTab$ensembl_gene_id %in% dummy,]
		resData[i,3]<-mean(dummy2$accEvents)
	}
	save(resData,file="data/resMean.Rdata")
}

# Pathway resampling table - Tab. S6
pathStats<-data.frame(pathway=pathways,meanPREs=double(length(pathways)),
	P_meanPREs=double(length(pathways)),relAcc=double(length(pathways)),
	P_relAcc=double(length(pathways)),accEvents=double(length(pathways)),
	P_accEvents=double(length(pathways)))
for(i in 1:length(pathways)){
	dummy<-geneTab[geneTab$ensembl_gene_id %in% lst$ensembl_gene_id[lst$pathway==pathways[i]],]
	pathStats$relAcc[i]<-mean(dummy$relAcc,na.rm=T)
	pathStats$P_relAcc[i]<-sum(resData$meanRelAcc>pathStats$relAcc[i])/n
	dummy<-genePREtab[genePREtab$ensembl_gene_id %in% lst$ensembl_gene_id[lst$pathway==pathways[i]],]
	pathStats$meanPREs[i]<-mean(dummy$PREs,na.rm=T)
	pathStats$P_meanPREs[i]<-sum(resData$meanNoPREs>pathStats$meanPREs[i])/n
	dummy<-uGeneTab[uGeneTab$ensembl_gene_id %in% lst$ensembl_gene_id[lst$pathway==pathways[i]],]
	pathStats$accEvents[i]<-mean(dummy$accEvents,na.rm=T)
	pathStats$P_accEvents[i]<-sum(resData$meanNoAcc>pathStats$accEvents[i])/n
}
tmp<-p.adjust(c(pathStats$P_meanPREs,pathStats$P_relAcc,pathStats$P_accEvents),method="BH")
pathStats$Q_meanPREs<-tmp[1:10]
pathStats$Q_relAcc<-tmp[11:20]
pathStats$Q_accEvents<-tmp[21:30]
write.table(pathStats,"pathwayStats.tsv",quote=F,sep="\t",row.names=F)

plot(resData$meanRelAcc,resData$meanNoPREs,ylim=c(5,45))
points(pathStats$relAcc,pathStats$meanPREs,col=2) # no

# Pathway resampling histogram plots - Fig. 4E,F
pdf("pathways-res.pdf"))
par(mfrow=c(3,1),las=1)
hist(resData$meanNoPREs,breaks=6:44,xlab="PREs per gene",ylab="Count",main="")
abline(v=pathStats$meanPREs)
abline(v=quantile(resData$meanNoPREs,.98),lwd=2,lty=2)
text(x=pathStats$meanPREs,y=1000,labels=pathStats$pathway,srt=45)
hist(resData$meanRelAcc,breaks=seq(.002,.012,.0001),
	xlab="Acceleration (%)",ylab="Count",main="")
abline(v=pathStats$relAcc)
abline(v=quantile(resData$meanRelAcc,.98),lwd=2,lty=2)
text(x=pathStats$relAcc,y=300,labels=pathStats$pathway,srt=45)
hist(resData$meanNoAcc,breaks=seq(2.5,23.5,.5),
	xlab="Acceleration events per gene",ylab="Count",main="")
abline(v=pathStats$accEvents)
abline(v=quantile(resData$meanNoAcc,.98),lwd=2,lty=2)
text(x=pathStats$accEvents,y=600,labels=pathStats$pathway,srt=45)
dev.off()

q(save="no")
