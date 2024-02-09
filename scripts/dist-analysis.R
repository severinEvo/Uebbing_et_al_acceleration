setwd("data/meme/")
median_na<-function(x){median(x,na.rm=T)}
quartile_na<-function(x){quantile(x,probs=c(.25,.75),na.rm=T)}
load("annotations/branch-list.Rdata")

allFiles<-list.files()

signFiles<-allFiles[grepl("fimo-sign-",allFiles)]
signTab<-read.delim(signFiles[1])
for(f in 2:length(signFiles)){
	tmp<-read.delim(signFiles[f])
	if(nrow(tmp)>0){
		signTab<-rbind(signTab,tmp)
}}
signTab<-signTab[!is.na(signTab$crossD),]

nsignFiles<-allFiles[grepl("fimo-nsign-",allFiles)]
nsignTab<-read.delim(nsignFiles[1])
for(f in 2:length(nsignFiles)){
	tmp<-read.delim(nsignFiles[f])
	if(nrow(tmp)>0){
		nsignTab<-rbind(nsignTab,tmp)
}}
nsignTab<-nsignTab[!is.na(nsignTab$crossD),]
tmpbra<-unique(nsignTab$branch)
tmpbra<-tmpbra[tmpbra %in% signTab$branch]
n<-length(tmpbra)
cDstats<-data.frame(branch=tmpbra,median_sign=double(n),median_nsign=double(n),
	wmu=double(n),wmu_P=double(n))
for(i in 1:n){
	tmpsign<-signTab$crossD[signTab$branch==cDstats[i,1]]
	tmpnsign<-nsignTab$crossD[nsignTab$branch==cDstats[i,1]]
	cDstats[i,2]<-median(tmpsign)
	cDstats[i,3]<-median_na(tmpnsign)
	tmpwmu<-wilcox.test(tmpsign,tmpnsign)
	cDstats[i,4]<-tmpwmu$statistic
	cDstats[i,5]<-tmpwmu$p.value
}
write.table(cDstats,"fimo-nsign-stats.tsv",sep="\t",quote=F,row.names=F)

signComp<-signTab[signTab$branch %in% cDstats$branch,]
nsignComp<-nsignTab[nsignTab$branch %in% cDstats$branch,]

# significance test for Fig. 3B
wilcox.test(nsignComp$crossD,signComp$crossD)
# W = 1.3491e+11, p-value = 0

pdf("../../acc-v-not.pdf") # Fig. 3B
plot(1,type="n",xlim=0:1,ylim=c(.5,2.5),bty="n",ylab="",yaxt="n",
	xlab="TF set dissimilarity between ingroup and outgroup")
axis(2,at=1:2,labels=c("No acceleration","Accelerated"),tick=F)
segments(x0=c(quantile(nsignComp$crossD,.05,na.rm=T),quantile(signComp$crossD,.05,na.rm=T)),y0=1:2,
	x1=c(quantile(nsignComp$crossD,.95,na.rm=T),quantile(signComp$crossD,.95,na.rm=T)),lwd=2)
rect(ytop=(1:2)-.1,xright=c(quartile_na(nsignComp$crossD)[2],quartile_na(signComp$crossD)[2]),
	ybottom=(1:2)+.1,xleft=c(quartile_na(nsignComp$crossD)[1],quartile_na(signComp$crossD)[1]),
	lwd=2,col=gray(.8))
segments(y0=(1:2)-.1,y1=(1:2)+.1,x0=c(median_na(nsignComp$crossD),median_na(signComp$crossD)),
	lwd=3)
dev.off()

q(save="no")
