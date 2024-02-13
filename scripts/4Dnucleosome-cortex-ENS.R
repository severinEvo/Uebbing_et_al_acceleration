reqPackages<-c("biomaRt")
newPackages<-reqPackages[!(reqPackages %in% utils::installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(newPackages)
library(biomaRt)

cortex<-read.delim("cortex_hiccups_loops_hg38.bed.gz",stringsAsFactors=F,header=F)
chroms<-sub("chr","",unique(cortex$V1))
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
allGenes<-getBM(filters="chromosome_name",
	attributes=c("ensembl_gene_id","external_gene_name","mmusculus_homolog_ensembl_gene"),
	values=list(chroms),mart=ensembl)
cortex<-merge(cortex,allGenes,by.x="V4",by.y="mmusculus_homolog_ensembl_gene")
write.table(cortex[,c(2:5)],gzfile("cortex_hiccups_loops_hg38_ENS.bed.gz"),quote=F,sep="\t",
	col.names=F,row.names=F)

q(save="no")
