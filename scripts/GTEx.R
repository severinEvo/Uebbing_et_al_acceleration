setwd("data/")
median_na<-function(x){median(x,na.rm=T)}
tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be set to 0.')
  if(any(x<0)) stop('Negative input values not permitted.')
  t<-sum(1-x/max(x))/(length(x)-1)
}

gtex<-read.delim("../annotations/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
gtex[3:ncol(gtex)]<-replace(gtex[3:ncol(gtex)],gtex[3:ncol(gtex)]==0,NA)
gtex[3:ncol(gtex)]<-log2(gtex[3:ncol(gtex)])
gtex[3:ncol(gtex)]<-apply(gtex[3:ncol(gtex)],2,z_fpkm)

transl<-read.delim("../annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
transl$SMTSD<-gsub(" - ","_",transl$SMTSD)
transl$SMTSD<-gsub(" ","-",transl$SMTSD)
tissues<-unique(transl$SMTSD)

median_gtex<-data.frame(matrix(ncol=length(tissues)+2,nrow=nrow(gtex)))
median_gtex[,1:2]<-gtex[,1:2]
colnames(median_gtex)<-c("Name","Description",tissues)
for(t in 1:length(tissues)){
	samples<-transl$SAMPID[transl$SMTSD==tissues[t]]
	samples<-gsub("-",".",samples)
	tmp<-gtex[,colnames(gtex) %in% samples]
	median_gtex[t+2]<-apply(tmp,1,median_na)
}
median_gtex[,3:ncol(median_gtex)]<-2^median_gtex[,3:ncol(median_gtex)]
median_gtex[,3:ncol(median_gtex)]<-replace(median_gtex[,3:ncol(median_gtex)],
	is.na(median_gtex[,3:ncol(median_gtex)]),0)

tau_gtex<-data.frame(matrix(ncol=3,nrow=nrow(gtex)))
tau_gtex[,1:2]<-gtex[,1:2]
colnames(tau_gtex)<-c("ensembl_gene_id","external_gene_name","tau")
tau_gtex$tau<-apply(median_gtex[,3:ncol(median_gtex)],1,tau)
tau_gtex$ensembl_gene_id<-substr(tau_gtex$ensembl_gene_id,1,15)
write.table(tau_gtex,"GTEx_tau.tsv",sep="\t",quote=F,row.names=F)

q(save="no")
