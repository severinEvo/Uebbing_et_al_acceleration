setwd("data/phyloP_out/")

dat<-NULL
fileList<-list.files()
phyPout<-fileList[grepl("[a-z|-]*_phyloP-phastBias.tsv.gz",fileList)]
for(i in phyPout){
	tmp<-read.delim(gzfile(i))[c(1:3,8,9)]
	dat<-rbind(dat,tmp)
}
pres<-read.delim("../PREs.bed.gz",header=F)
colnames(pres)<-c("chr","start","end","PRE_type","seq_type","phastCons")
datn<-merge(dat,pres)
qual<-read.delim("../qualFilter/qualFilter_useBranches.bed.gz")
colnames(qual2)<-colnames(qual)
qual<-rbind(qual,qual2)
datq<-merge(datn,qual)
datq<-unique(datq)
write.table(datq,gzfile("phyloP-all.tsv.gz"),quote=F,sep="\t",row.names=F)

sign<-datq[p.adjust(datq$pval,method="BH")<.05,]
write.table(sign[,1:5],gzfile("phyloP-sign.bed.gz"),sep="\t",row.names=F,col.names=F,quote=F)
write.table(sign[,1:5],gzfile("phyloP-sign_0based.tsv.gz"),sep="\t",row.names=F,quote=F)

# Basic PRE stats
usign<-unique(sign[,c(1:3,6:8)])
udatq<-unique(datq[,c(1:3,6:8)])
nrow(udatq) # 434,240 - Total # of usable PREs
# 319,227 pCEs / 115,013 cCREs
# 118,008 intergenic / 84,082 promoter-proximal / 232,150 intronic
write.table(udatq,gzfile("../PRE-gene-dataset/PREs_qualFiltered.bed.gz"),
	col.names=F,row.names=F,quote=F,sep="\t")
table(usign$PRE_type)/table(udatq$PRE_type)
#      cCRE       pCE
# 0.3117039 0.1420463
table(sign$branch)/table(datq$branch)

ccres<-udatq[udatq$PRE_type=="cCRE",]
pces<-udatq[udatq$PRE_type=="pCE",]
write.table(pces[,1:3],gzfile("../PRE-gene-dataset/pCEs-qualFiltered.bed.gz"),
	quote=F,sep="\t",col.names=F,row.names=F)
write.table(ccres[,1:3],gzfile("../PRE-gene-dataset/cCREs-qualFiltered.bed.gz"),
	quote=F,sep="\t",col.names=F,row.names=F)

load("../../annotations/branch-list.Rdata")
branches<-names(branchlist)
branches<-branches[branches %in% unique(sign$branch)]
signTab<-data.frame(table(sign$branch))
colnames(signTab)<-c("branch","sign")
signTab$branch<-as.character(signTab$branch)
testTab<-data.frame(table(datq$branch))
colnames(testTab)<-c("branch","tested")
testTab$branch<-as.character(testTab$branch)
tab<-merge(signTab,testTab)
tab$relAcc<-tab$sign/tab$tested
o<-ordered(tab$branch,levels=branches)
pdf("phyP-sign-hist.pdf",width=20)
barplot(tab$sign ~ o,las=3,xlab="",ylab="Accelerated PREs")
dev.off()
pdf("phyP-relAcc-hist.pdf",width=20)
barplot(tab$relAcc ~ o,las=3,xlab="",ylab="Accelerated PREs (%)")
dev.off()

branchlen<-read.delim("../../annotations/branch-lens.tsv",)
tab<-merge(tab,branchlen)
cols<-c(rep("#b2df8a",20),0,0,rep("#33a02c",16),0,0,0,rep("#a6cee3",2),rep("#1f78b4",6),
	rep("#a6cee3",11),0,0,0,rep("#e31a1c",12),0,0,rep("#fb9a99",7),rep(0,5),rep("#fdbf6f",6),
	0,0,0,rep("#ff7f00",2))

cor.test(tab$relAcc,tab$branch_len,method="spearman")
# S = 91844, p-value = 3.616e-06, rho = 0.4488809

# Branch length vs acceleration frequency (Fig 1C)
pdf("../../branchlength-v-freqAcc.pdf")
plot(tab$relAcc*100,tab$branch_len,las=1,bg=cols,pch=21,cex=1.5,bty="n",
	xlab="Accelerated PREs (%)",ylab="Branch length (substitutions/site)")
legend("topleft",c(expression(paste(rho," = 0.45")),expression(paste("P = 3.6 ",10^-6))),bty="n")
plot(tab$relAcc*100,tab$branch_len,type="n",
	xlab="Accelerated PREs (%)",ylab="Branch length (substitutions/site)")
text(tab$relAcc*100,tab$branch_len,tab$branch)
dev.off()

# Relative acceleratio per branch and PRE type (Fig 1B)
signTab<-data.frame(table(sign[,c(4,6)]))
colnames(signTab)[3]<-"sign"
signTab$branch<-as.character(signTab$branch)
signTab$PRE_type<-as.character(signTab$PRE_type)
testTab<-data.frame(table(datq[,c(4,6)]))
colnames(testTab)[3]<-"tested"
testTab$branch<-as.character(testTab$branch)
testTab$PRE_type<-as.character(testTab$PRE_type)
tab<-merge(signTab,testTab)
tab$relAcc<-tab$sign/tab$tested
tab$name<-paste(tab$branch,tab$PRE_type,sep="_")
o<-ordered(tab$branch[tab$PRE_type=="pCE"],levels=branches)
pdf("../../phyP-mirror-PRE-relAcc-hist.pdf",width=20)
barplot(tab$relAcc[tab$PRE_type=="pCE"] ~ o,las=3,ylim=c(0,.06),yaxt="n",
	xlab="",ylab="Accelerated pCEs (%)")
axis(2,(0:6)/100,0:6)
barplot(tab$relAcc[tab$PRE_type=="cCRE"] ~ o,las=3,ylim=c(0,.06),yaxt="n",
	xlab="",ylab="Accelerated cCREs (%)")
axis(2,(0:6)/100,0:6)
dev.off()

# Filtering effects on PREs (Fig S2)
datn$pos<-paste(datn$chr,":",datn$start,"-",datn$end,sep="")
datq$pos<-paste(datq$chr,":",datq$start,"-",datq$end,sep="")
sign<-datq[p.adjust(datq$pval,method="BH")<.05,]
pdf("../../PRE-type-filter-dist.pdf")
par(mfrow=c(3,2),las=1)
hist(table(datn$pos[datn$PRE_type=="pCE"]),main="Unique pCEs",xlab="",xaxt="n",ylab="",
	breaks=0:100)
hist(table(datn$pos[datn$PRE_type=="cCRE"]),main="Unique cCREs",xlab="",xaxt="n",ylab="",
	breaks=0:100)
hist(table(datq$pos[datq$PRE_type=="pCE"]),main="QualityFiltered pCEs",xlab="",xaxt="n",
	ylab="Number of PREs",breaks=0:100)
hist(table(datq$pos[datq$PRE_type=="cCRE"]),main="QualityFiltered cCREs",xlab="",xaxt="n",
	ylab="",breaks=0:100)
hist(table(sign$pos[sign$PRE_type=="pCE"]),main="Accelerated pCEs",xlab="Number of branches",
	ylab="",breaks=0:100)
hist(table(sign$pos[sign$PRE_type=="cCRE"]),main="Accelerated cCREs",xlab="Number of branches",
	ylab="",breaks=0:100)
dev.off()
pdf("../../seq-type-filter-dist.pdf")
par(mfrow=c(3,3),las=1)
hist(table(datn$pos[datn$seq_type=="promoter-proximal"]),main="promoter-proximal",xlab="",
	xaxt="n",ylab="",breaks=0:100)
hist(table(datn$pos[datn$seq_type=="intronic"]),main="intronic",xlab="",
	xaxt="n",ylab="",breaks=0:100)
hist(table(datn$pos[datn$seq_type=="intergenic"]),main="intergenic",xlab="",
	xaxt="n",ylab="",breaks=0:100)
hist(table(datq$pos[datq$seq_type=="promoter-proximal"]),main="",xlab="",
	xaxt="n",ylab="Number of PREs",breaks=0:100)
hist(table(datq$pos[datq$seq_type=="intronic"]),main="QualityFiltered",xlab="",
	xaxt="n",ylab="",breaks=0:100)
hist(table(datq$pos[datq$seq_type=="intergenic"]),main="",xlab="",
	xaxt="n",ylab="",breaks=0:100)
hist(table(sign$pos[sign$seq_type=="promoter-proximal"]),main="",xlab="",
	xaxt="n",ylab="",breaks=0:100)
hist(table(sign$pos[sign$seq_type=="intronic"]),main="Accelerated",xlab="Number of branches",
	xaxt="n",ylab="",breaks=0:100)
hist(table(sign$pos[sign$seq_type=="intergenic"]),main="",xlab="",
	xaxt="n",ylab="",breaks=0:100)
dev.off()

# BED file for UCSC genome browser:
# chr1    937873  937930  ACC:pCE:tested  100     .
# chr1    966648  966798  ACC:cCRE:afrosoricida   1000    .
beds<-merge(udatq[,1:4],sign[,1:4],all.x=T)
beds$name<-character(nrow(beds))
beds$name[is.na(beds$branch)]<-paste(beds$PRE_type[is.na(beds$branch)],"tested",sep=":")
beds$name[!is.na(beds$branch)]<-paste("ACC",beds$PRE_type[!is.na(beds$branch)],
	beds$branch[!is.na(beds$branch)],sep=":")
beds$score<-1000
beds$score[is.na(beds$branch)]<-100
beds$strand<-"."
write.table(beds[beds$chr!="chrY",-c(4:5)],"acceleration.bed",
	col.names=F,row.names=F,quote=F,sep="\t")

q(save="no")
