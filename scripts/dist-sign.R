#!/usr/bin/env Rscript
setwd("data/meme/")
"%out%"<-function(x,table){match(x,table,nomatch=0)==0}
all<-read.delim("../../phyloP_out/phyloP-all.tsv.gz")
sign<-all[p.adjust(all$pval,method="BH")<.05,]
load("../../../annotations/branch-list.Rdata")

args<-commandArgs(trailingOnly=T)
branch<-args[1]
rm(all)
sign<-sign[sign$branch==branch,]

allfiles<-list.files()
usefiles<-allfiles[grepl("fimo_all_chr.*.bed",allfiles)]

sign$inD<-NA
sign$outD<-NA
sign$allD<-NA
sign$crossD<-NA
bpos<-which(names(branchlist)==branch)
brdummy<-data.frame(branch=branch,spec=unlist(clalist[[bpos]]),
	group=c(rep("in",length(unlist(clalist[[bpos]][1:2]))),
		rep("out",length(unlist(clalist[[bpos]][3])))))
brSign<-sign[,1:5]
brSign$filename<-paste(brSign$chr,"/",brSign$chr,":",brSign$start,"-",brSign$end,"_fimo_all_",brSign$chr,".bed",sep="")
spSign<-merge(sign[,1:5],brdummy)
spSign$filename<-paste(spSign$chr,"/",spSign$chr,":",spSign$start,"-",spSign$end,"_fimo_all_",spSign$chr,".bed",sep="")
for(i in 1:nrow(brSign)){
	if(brSign$filename[i] %in% usefiles){
		inbed<-read.delim(brSign$filename[i],stringsAsFactors=F,header=F)
	}else{
		next
	}
	tspSign<-spSign[spSign$start==brSign$start[i] & spSign$end==brSign$end[i] & spSign$spec %in% inbed$V7,]
	in_tspSign<-tspSign[tspSign$group=="in",]
	out_tspSign<-tspSign[tspSign$group=="out",]
	if(nrow(in_tspSign)>1 & nrow(out_tspSign)>1){
		for(j in 1:nrow(tspSign)){
			tmp<-inbed[inbed$V7==tspSign$spec[j],]
			if(nrow(tmp)>0){
				assign(tmp$V7[1],data.frame(table(tmp$V4)))
		}}
		in_table<-merge(get(in_tspSign$spec[1]),get(in_tspSign$spec[2]),by="Var1",all=T)
		colnames(in_table)[1:3]<-c("TF",in_tspSign$spec[1:2])
		if(nrow(in_tspSign)>2){
			for(k in 3:nrow(in_tspSign)){
				in_table<-merge(in_table,get(in_tspSign$spec[k]),by.x="TF",by.y="Var1",all=T)
				colnames(in_table)[k+1]<-in_tspSign$spec[k]
		}}
		in_table[,2:ncol(in_table)]<-replace(in_table[,2:ncol(in_table)],
			is.na(in_table[,2:ncol(in_table)]),0)
		in_d<-dist(t(in_table[,2:ncol(in_table)]),method="binary")
		out_table<-merge(get(out_tspSign$spec[1]),get(out_tspSign$spec[2]),by="Var1",all=T)
		colnames(out_table)[1:3]<-c("TF",out_tspSign$spec[1:2])
		if(nrow(out_tspSign)>2){
			for(k in 3:nrow(out_tspSign)){
				out_table<-merge(out_table,get(out_tspSign$spec[k]),by.x="TF",by.y="Var1",all=T)
				colnames(out_table)[k+1]<-out_tspSign$spec[k]
		}}
		out_table[,2:ncol(out_table)]<-replace(out_table[,2:ncol(out_table)],
			is.na(out_table[,2:ncol(out_table)]),0)
		out_d<-dist(t(out_table[,2:ncol(out_table)]),method="binary")
		table<-merge(in_table,out_table,all=T)
		table[,2:ncol(table)]<-replace(table[,2:ncol(table)],
			is.na(table[,2:ncol(table)]),0)
		all_d<-dist(t(table[,2:ncol(table)]),method="binary")
		sign$inD[sign$start==tspSign$start[1] & sign$end==tspSign$end[1]]<-mean(in_d)
		sign$outD[sign$start==tspSign$start[1] & sign$end==tspSign$end[1]]<-mean(out_d)
		sign$allD[sign$start==tspSign$start[1] & sign$end==tspSign$end[1]]<-mean(all_d)
		sign$crossD[sign$start==tspSign$start[1] & sign$end==tspSign$end[1]]<-mean(as.vector(all_d)[as.vector(all_d) %out% c(as.vector(in_d),as.vector(out_d))])
}}
write.table(sign[!is.na(sign$inD),],
	gzfile(paste("fimo-sign-",branch,".tsv.gz",sep="")),
	sep="\t",quote=F,row.names=F)

q(save="no")
