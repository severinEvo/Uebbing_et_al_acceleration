quartile<-function(x){quantile(x,probs=c(.25,.75))}
quantile90<-function(x){quantile(x,probs=c(.05,.95))}

pres<-read.delim("data/PRE-gene-dataset/PREs_qualFiltered.bed.gz",header=F)
colnames(pres)<-c("chr","start","end","PRE_type","seq_type","phastCons")
pres$pos<-paste(pres$chr,":",pres$start,"-",pres$end,sep="")

HiC<-read.delim("data/PRE-gene-dataset/connection-collection.bed.gz",header=F)
colnames(HiC)<-c("chr","start","end","ensembl_gene_id","PRE_type","connection")
HiC$pos<-paste(HiC$chr,":",HiC$start,"-",HiC$end,sep="")
uHiC<-unique(HiC[,c("ensembl_gene_id","PRE_type","pos")])
if(file.exists("data/allGenesGO.tsv.gz")){
	allGenesGO<-read.delim("data/allGenesGO.tsv.gz")
}else{
	chroms<-sub("chr","",unique(uHiC$chr))
	library(biomaRt)
	ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
	allGenesGO<-getBM(
		attributes=c("ensembl_gene_id","external_gene_name","go_id","name_1006","namespace_1003"),
		filters=c("ensembl_gene_id","chromosome_name","with_go"),
		values=list(unique(uHiC$ensembl_gene_id),chroms,T),mart=ensembl)
	write.table(allGenesGO,gzfile("data/allGenesGO.tsv.gz"),quote=F,sep="\t",row.names=F)
}
goodGO<-table(allGenesGO$go_id)
goodGO<-goodGO[goodGO >5]
allGenesGO<-allGenesGO[allGenesGO$go_id %in% names(goodGO),]

datq<-read.delim("data/phyloP_out/phyloP-all.tsv.gz")
datq$pos<-paste(datq$chr,":",datq$start,"-",datq$end,sep="")
testTab<-data.frame(table(datq$pos))
colnames(testTab)<-c("pos","tested")
testTab$pos<-as.character(testTab$pos)
sign<-datq[p.adjust(datq$pval,method="BH")<.05,]
sign$sign<-1
signTab<-data.frame(table(sign$pos))
colnames(signTab)<-c("pos","sign")
signTab$pos<-as.character(signTab$pos)
signTab<-merge(testTab,signTab,all.x=T)
signTab$sign<-replace(signTab$sign,is.na(signTab$sign),0)
geneData<-merge(unique(uHiC[,c("ensembl_gene_id","pos")]),signTab)

# cCREs vs pCEs - Tab. S1
tmp<-merge(uHiC[,c("ensembl_gene_id","pos")],allGenesGO[,c("ensembl_gene_id","go_id")])
bg<-data.frame(table(tmp$go_id))
colnames(bg)<-c("go_id","bg")
tmp<-merge(uHiC[uHiC$PRE_type=="pCE",c("ensembl_gene_id","pos")],
	allGenesGO[,c("ensembl_gene_id","go_id")])
pces<-data.frame(table(tmp$go_id))
colnames(pces)<-c("go_id","pCEs")
tmp<-merge(uHiC[uHiC$PRE_type=="cCRE",c("ensembl_gene_id","pos")],
	allGenesGO[,c("ensembl_gene_id","go_id")])
ccres<-data.frame(table(tmp$go_id))
colnames(ccres)<-c("go_id","cCREs")
all<-merge(bg,pces,all.x=T)
all<-merge(all,ccres,all.x=T)

q<-all$pCEs # number of white balls drawn # sign of GO x
m<-sum(all$pCEs) # number of white balls in the urn. # all sign acc
n<-sum(all$bg)-sum(all$pCEs) # number of black balls in the urn. # all non-sign (non-acc)
k<-all$bg # number of balls drawn from the urn. # bg of GO x
all$P_pCE<-phyper(q-1,m,n,k,lower.tail=F)
q<-all$cCREs # number of white balls drawn # sign of GO x
m<-sum(all$cCREs) # number of white balls in the urn. # all sign acc
n<-sum(all$bg)-sum(all$cCREs) # number of black balls in the urn. # all non-sign (non-acc)
all$P_cCRE<-phyper(q-1,m,n,k,lower.tail=F)
all$enrichment_pCE<-all$pCEs/sum(all$pCEs)/all$bg*sum(all$bg)
all$enrichment_cCRE<-all$cCREs/sum(all$cCREs)/all$bg*sum(all$bg)
resTab<-merge(unique(allGenesGO[,3:5]),all[,c(1,5:8)],all.x=T)
pTab<-data.frame(pivot_longer(resTab,cols=4:5,names_to="PRE_type",values_to="P"))
pTab$PRE_type<-sub("P_","",pTab$PRE_type)
pTab$Q<-p.adjust(pTab$P,method="BH")
pTabsign<-pTab[pTab$Q<.05 & !is.na(pTab$P),]
write.table(pTabsign,"PRE-GO-sign.tsv",sep="\t",quote=F,row.names=F)

cntTab<-data.frame(table(uHiC$ensembl_gene_id))
colnames(cntTab)<-c("ensembl_gene_id","connPREs")
cntTab$ensembl_gene_id<-as.character(cntTab$ensembl_gene_id)
cntTab<-merge(cntTab,unique(allGenesGO[,1:2]))
cntTab<-data.frame(cntTab,tested=double(nrow(cntTab)),sign=double(nrow(cntTab)))
for(i in 1:nrow(cntTab)){
	tmp<-geneData[geneData$ensembl_gene_id==cntTab$ensembl_gene_id[i],]
	cntTab[i,4:5]<-colSums(tmp[,3:4])
}
cntTab<-cntTab[cntTab$tested>0,]
cntTab$pcAcc<-cntTab$sign/cntTab$tested*100
htfs<-read.delim("annotations/human-tfs_1-s2-0-S0092867418301065-mmc2_edit.tsv")

# Mann-Whitney U tests for TF gene against non-TF gene stats (Fig. 4B-D)
wilcox.test(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F],
	cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])$p.value
# 1.871107e-44
summary(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])
summary(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])
wilcox.test(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F],
	cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])$p.value
# 4.819392e-05
summary(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])
summary(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])
wilcox.test(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F],
	cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])$p.value
# 1.215249e-44
summary(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])
summary(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])

# Plot Fig. 4B-D
pdf("TFs-v-other.pdf",height=4)
par(mfrow=c(1,3),las=1)
plot(1,type="n",xlim=c(-.4,1.4),ylim=c(0,150),xlab="",xaxt="n",ylab="PREs per gene",frame.plot=F)
segments(x0=0:1,y0=c(
	quantile90(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
	quantile90(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1]),y1=c(
	quantile90(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
	quantile90(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
segments(x0=rep((0:1)-.2,2),x1=rep((0:1)+.2,2),y0=c(
	quantile90(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
	quantile90(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1],
	quantile90(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
	quantile90(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
rect(col=gray(.8),xleft=(0:1)-.4,xright=(0:1)+.4,ybottom=c(
		quartile(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
		quartile(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1]),ytop=c(
		quartile(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
		quartile(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
segments(x0=(0:1)-.4,x1=(0:1)+.4,y0=c(
		median(cntTab$connPREs[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]]),
		median(cntTab$connPREs[(cntTab$ensembl_gene_id %in% htfs$ID)==F])),lwd=2)
text(.5,140,expression(paste("P = 1.9"^"-44")))
axis(1,at=0:1,labels=c("Transcription factors","Other genes"),tick=F)

plot(1,type="n",xlim=c(-.4,1.4),ylim=c(0,2.2),xlab="",xaxt="n",ylab="Acceleration (%)",
	frame.plot=F)
segments(x0=0:1,y0=c(
	quantile90(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
	quantile90(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1]),y1=c(
	quantile90(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
	quantile90(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
segments(x0=rep((0:1)-.2,2),x1=rep((0:1)+.2,2),y0=c(
	quantile90(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
	quantile90(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1],
	quantile90(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
	quantile90(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
rect(col=gray(.8),xleft=(0:1)-.4,xright=(0:1)+.4,ybottom=c(
		quartile(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
		quartile(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1]),ytop=c(
		quartile(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
		quartile(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
segments(x0=(0:1)-.4,x1=(0:1)+.4,y0=c(
		median(cntTab$pcAcc[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]]),
		median(cntTab$pcAcc[(cntTab$ensembl_gene_id %in% htfs$ID)==F])),lwd=2)
text(.5,2.0,expression(paste("P = 4.8"^"-5")))
axis(1,at=0:1,labels=c("Transcription factors","Other genes"),tick=F)

plot(1,type="n",xlim=c(-.4,1.4),ylim=c(0,55),xlab="",xaxt="n",
	ylab="Acceleration events per gene",frame.plot=F)
segments(x0=0:1,y0=c(
	quantile90(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
	quantile90(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1]),y1=c(
	quantile90(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
	quantile90(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
segments(x0=rep((0:1)-.2,2),x1=rep((0:1)+.2,2),y0=c(
	quantile90(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
	quantile90(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1],
	quantile90(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
	quantile90(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
rect(col=gray(.8),xleft=(0:1)-.4,xright=(0:1)+.4,ybottom=c(
		quartile(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[1],
		quartile(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[1]),ytop=c(
		quartile(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]])[2],
		quartile(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])[2]))
segments(x0=(0:1)-.4,x1=(0:1)+.4,y0=c(
		median(cntTab$sign[cntTab$ensembl_gene_id %in% htfs$ID[htfs$Is.it.TF.=="Yes"]]),
		median(cntTab$sign[(cntTab$ensembl_gene_id %in% htfs$ID)==F])),lwd=2)
text(.5,50,expression(paste("P = 1.2"^"-44")))
axis(1,at=0:1,labels=c("Transcription factors","Other genes"),tick=F)
dev.off()

# Branch-specific (percent acc) - Tab. S2
load("annotations/branch-list.Rdata")
branches<-names(branchlist)

bigTab<-merge(allGenesGO[,c(1,3)],unique(HiC[,c(4,7)]))
resT<-unique(allGenesGO[,3:4])
for(b in 1:length(branches)){
	print(paste(b,"/",length(branches),"branches"))
	tmp2<-merge(sign[sign$branch==branches[b],6:7],bigTab,all.y=T)
	tmp2$sign<-replace(tmp2$sign,is.na(tmp2$sign),0)
	bg<-data.frame(table(tmp2$go_id))
	colnames(bg)<-c("go_id","bg")
	bg$go_id<-as.character(bg$go_id)
	if(sum(tmp2$sign)==0){
		all<-data.frame(bg,acc=0)
	}else{
		acc<-data.frame(table(tmp2$go_id[tmp2$sign==1]))
		colnames(acc)<-c("go_id","acc")
		acc$go_id<-as.character(acc$go_id)
		all<-merge(bg,acc,all.x=T)
		all$acc<-replace(all$acc,is.na(all$acc),0)
	}
	q<-all$acc # number of white balls drawn # sign of GO x
	m<-sum(all$acc) # number of white balls in the urn. # all sign acc
	n<-sum(all$bg)-sum(all$acc) # number of black balls in the urn. # all non-acc
	k<-all$bg # number of balls drawn from the urn. # bg of GO x
	all$P<-phyper(q-1,m,n,k,lower.tail=F)
	colnames(all)[4]<-paste("P",branches[b],sep="_")
	resT<-merge(resT,all[,c(1,4)],all.x=T)
}
pTab<-data.frame(pivot_longer(resT,cols=3:102,names_to="branch",values_to="P"))
pTab$branch<-sub("P_","",pTab$branch)
pTab$Q<-p.adjust(pTab$P,method="BH")
pTabsign<-pTab[pTab$Q<.05 & !is.na(pTab$P),]
write.table(pTabsign,"all_GObiol-hyperg-sign.tsv",quote=F,row.names=F,sep="\t")

# Figure S4
GOslimagr<-read.delim("annotations/GOslim_agr.tsv") # Downloaded from QuickGO
GOslim<-unique(GOslimagr[,5:6])
colnames(GOslim)<-c("goslim_id","go_id")
goreslim<-merge(GOslim,pTabsign,all.y=T)
goreslim<-rbind(goreslim,data.frame(pTabsign,goslim_id="GO:0008150")[,c(1,6,2:4)])

branchlist<-names(branchlist)[names(branchlist) %in% unique(goreslim$branch)]
goslimsort<-read.delim("annotations/goslimlist-sorted.tsv") # Downloaded from QuickGO

library(RColorBrewer)
palette<-brewer.pal(9,'Blues')
pdf("GOslim-palette.pdf",height=25)
par(las=2)
plot(1,type="n",xlim=c(2,22),ylim=c(3,83),xlab="",xaxt="n",ylab="",yaxt="n",bty="n")
axis(2,at=85:1,labels=branchlist,tick=F)
axis(3,at=2:22,labels=goslimsort$goslim_description[-1],tick=F)
for(b in 1:length(branchlist)){
	branch<-branchlist[b]
	tmptab<-data.frame(table(goreslim$goslim_id[goreslim$branch==branch]))
	colnames(tmptab)<-c("goslim_id","count")
	tmptab<-merge(goslimsort,tmptab,all=T)
	tmptab$count<-replace(tmptab$count,is.na(tmptab$count),0)
	tmptab<-tmptab[match(goslimsort$goslim_id,tmptab$goslim_id),]
	cols<-palette[floor(log2(tmptab$count+1)/max(log2(tmptab$count+1))*8)+1]
	rect(xleft=c(0,2:22)-.4,ybottom=86-b-.4,xright=c(0,2:22)+.4,ytop=86-b+.4,col=cols)
}
dev.off()

# hyperG GO test for % accelerated PREs per genes - Tab. S3
tmp<-allGenesGO[allGenesGO$ensembl_gene_id %in% cntTab$ensembl_gene_id,c(1,3)]
bg<-data.frame(table(tmp$go_id))
colnames(bg)<-c("go_id","bg")
bg$go_id<-as.character(bg$go_id)
tmp<-allGenesGO[allGenesGO$ensembl_gene_id %in% 
	cntTab$ensembl_gene_id[cntTab$pcAcc >quantile(cntTab$pcAcc,.95)],c(1,3)]
top<-data.frame(table(tmp$go_id))
colnames(top)<-c("go_id","top")
top$go_id<-as.character(top$go_id)
all<-merge(bg,top,all.x=T)
all$top<-replace(all$top,is.na(all$top),0)
q<-all$top # number of white balls drawn # sign of GO x
m<-sum(all$top) # number of white balls in the urn. # all sign top
n<-sum(all$bg)-sum(all$top) # number of black balls in the urn. # all non-top
k<-all$bg # number of balls drawn from the urn. # bg of GO x
all$P<-phyper(q-1,m,n,k,lower.tail=F)
all$Q<-p.adjust(all$P,method="BH")
all$enrichment<-all$top/sum(all$top)/all$bg*sum(all$bg)
signTopGO<-merge(unique(allGenesGO[,-(1:2)]),all[all$Q<.05,])
signTopCats<-merge(tmp,allGenesGO[allGenesGO$go_id %in% signTopGO$go_id,])
write.table(signTopGO,"pcAcc-hyperG-sign.tsv",quote=F,sep="\t",row.names=F)

pdf("pcAcc-hist.pdf",width=10) # Like Fig. 4A but for % acceleration
hist(cntTab$pcAcc,seq(0,20,.1),las=1,xlim=c(0,10),
	xlab="Acceleration per PRE (%)",ylab="Count",main="")
abline(v=quantile(cntTab$pcAcc,.95))
dev.off()

# hyperG GO test for no. acc events - Tab. S5
tmp<-allGenesGO[allGenesGO$ensembl_gene_id %in% cntTab$ensembl_gene_id,c(1,3)]
bg<-data.frame(table(tmp$go_id))
colnames(bg)<-c("go_id","bg")
bg$go_id<-as.character(bg$go_id)
tmp<-allGenesGO[allGenesGO$ensembl_gene_id %in% 
	cntTab$ensembl_gene_id[cntTab$sign >quantile(cntTab$sign,.95)],c(1,3)]
top<-data.frame(table(tmp$go_id))
colnames(top)<-c("go_id","top")
top$go_id<-as.character(top$go_id)
all<-merge(bg,top,all.x=T)
all$top<-replace(all$top,is.na(all$top),0)
q<-all$top # number of white balls drawn # sign of GO x
m<-sum(all$top) # number of white balls in the urn. # all sign top
n<-sum(all$bg)-sum(all$top) # number of black balls in the urn. # all non-top
k<-all$bg # number of balls drawn from the urn. # bg of GO x
all$P<-phyper(q-1,m,n,k,lower.tail=F)
all$Q<-p.adjust(all$P,method="BH")
all$enrichment<-all$top/sum(all$top)/all$bg*sum(all$bg)
signTopGO<-merge(unique(allGenesGO[,-(1:2)]),all[all$Q<.05,])
signTopCats<-merge(tmp,allGenesGO[allGenesGO$go_id %in% signTopGO$go_id,])
write.table(signTopGO,"../noAccEventsPerGene-hyperG-sign.tsv",quote=F,sep="\t",row.names=F)

# Fig. 4A
pdf("noAccEvents.pdf",width=10)
hist(cntTab$sign,seq(0,265,1),las=1,xlim=c(0,60),ylim=c(0,1700),
	xlab="Acceleration events per gene",ylab="Count",main="")
abline(v=quantile(cntTab$sign,.95))
dev.off()

# hyperG GO test for no. connected PREs - Tab. S4
tmp<-allGenesGO[allGenesGO$ensembl_gene_id %in% cntTab$ensembl_gene_id,c(1,3)]
bg<-data.frame(table(tmp$go_id))
colnames(bg)<-c("go_id","bg")
bg$go_id<-as.character(bg$go_id)
tmp<-allGenesGO[allGenesGO$ensembl_gene_id %in% 
	cntTab$ensembl_gene_id[cntTab$connPREs >quantile(cntTab$connPREs,.95)],c(1,3)]
top<-data.frame(table(tmp$go_id))
colnames(top)<-c("go_id","top")
top$go_id<-as.character(top$go_id)
all<-merge(bg,top,all.x=T)
all$top<-replace(all$top,is.na(all$top),0)
q<-all$top # number of white balls drawn # sign of GO x
m<-sum(all$top) # number of white balls in the urn. # all sign top
n<-sum(all$bg)-sum(all$top) # number of black balls in the urn. # all non-top
k<-all$bg # number of balls drawn from the urn. # bg of GO x
all$P<-phyper(q-1,m,n,k,lower.tail=F)
all$Q<-p.adjust(all$P,method="BH")
all$enrichment<-all$top/sum(all$top)/all$bg*sum(all$bg)
signTopGO<-merge(unique(allGenesGO[,-(1:2)]),all[all$Q<.05,])
signTopCats<-merge(tmp,allGenesGO[allGenesGO$go_id %in% signTopGO$go_id,])
write.table(signTopGO,"noPREsPerGene-hyperG-sign.tsv",quote=F,sep="\t",row.names=F)

pdf("noPREsPerGene.pdf",width=10) # Like Fig. 4A but for #PREs per gene
hist(cntTab$connPREs,seq(0,503,1),las=1,xlim=c(0,200),
	xlab="PREs per gene",ylab="Count",main="")
abline(v=quantile(cntTab$connPREs,.95))
dev.off()

q(save="no")
