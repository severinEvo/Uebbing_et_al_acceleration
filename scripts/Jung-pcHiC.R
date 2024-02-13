#!/usr/bin/Rscript
reqPackages<-c("biomaRt")
newPackages<-reqPackages[!(reqPackages %in% utils::installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(newPackages)
library(biomaRt)

dat<-unique(read.delim("NPC_P-O-contacts_hg38_1perline.bed.gz",stringsAsFactors=F,header=F))
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
exprGenes<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),
	filters="external_gene_name",values=unique(dat$V4),mart=ensembl)
datb<-merge(dat,exprGenes,by.x="V4",by.y="external_gene_name")
write.table(datb[,2:5],gzfile("NPC_P-O-contacts_hg38_ENSEMBL.bed.gz"),quote=F,sep="\t",
	row.names=F,col.names=F)

q(save="no")
