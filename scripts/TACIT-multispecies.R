setwd("data/TACIT/")
quartile_na<-function(x){quantile(x,probs=c(.25,.75),na.rm=T)}
quantile95_na<-function(x){quantile(x,probs=c(.025,.975),na.rm=T)}
median_na<-function(x){median(x,na.rm=T)}

datq<-read.delim("../phyloP_out/phyloP-all.tsv.gz")
datq$pos<-paste(datq$chr,":",datq$start,"-",datq$end,sep="")
nsign<-datq[p.adjust(datq$pval,method="BH")>.05,]
sign<-datq[p.adjust(datq$pval,method="BH")<.05,]
load("../../annotations/zoonomia-branch-list.Rdata")
branches<-unique(sign$branch)
branches<-branches[branches %in% names(branchlist)]

cortex<-read.delim("cortex_ocrs_filtered_named_hg38_CRE-intersect.bed")
cortex$pos<-paste(cortex$chr,":",cortex$start,"-",cortex$end,sep="")
colnames(cortex)<-gsub(".","_",colnames(cortex),fixed=T)
cortexO<-read.delim("cortex_ocrs_filtered_named.txt.gz")
colnames(cortexO)<-gsub(".","_",colnames(cortexO),fixed=T)
cortex<-merge(cortex,cortexO)
cortex[6:ncol(cortex)]<-replace(cortex[6:ncol(cortex)],cortex[6:ncol(cortex)]==-1,NA)

cortex_tab<-data.frame(branch=branches,P=double(length(branches)),direction=logical(length(branches)),
	n_sign=double(length(branches)))
pdf("TACIT-cortex.pdf",height=4)
for(b in branches){
	bpos<-which(names(branchlist)==b)
	cortex_tab[bpos,1]<-b
	brdummy<-data.frame(branch=b,spec=unlist(branchlist[[bpos]]),
		group=c(rep("in",length(unlist(branchlist[[bpos]][1:2]))),
			rep("out",length(unlist(branchlist[[bpos]][3])))))
	brdummy<-brdummy[brdummy$spec %in% colnames(cortex),]
	if(nrow(brdummy[brdummy$group=="in",])<1 | nrow(brdummy[brdummy$group=="out",])<1){
		next
	}
	brSign<-sign[sign$branch==b,4:6]
	brNsign<-nsign[nsign$branch==b,c(4,5,9)]
	tmpSign<-cortex[cortex$pos %in% brSign$pos,]
	in_sign<-tmpSign[,brdummy$spec[brdummy$group=="in"]]
	cortex_tab[bpos,4]<-nrow(in_sign)
	if(is.data.frame(in_sign) && dim(in_sign)[2]>1){
		mean_inS<-rowMeans(in_sign,na.rm=T)
	}else if(is.vector(in_sign) && length(in_sign)>1){
		mean_inS<-in_sign
	}else{next}
	if(sum(mean_inS,na.rm=T)==0){next}
	out_sign<-tmpSign[,brdummy$spec[brdummy$group=="out"]]
	if(is.data.frame(out_sign) && dim(out_sign)[2]>1){
		mean_outS<-rowMeans(out_sign,na.rm=T)
	}else if(is.vector(out_sign) && length(out_sign)>1){
		mean_outS<-out_sign
	}else{next}
	if(sum(mean_outS,na.rm=T)==0){next}
	tmpNsign<-cortex[cortex$pos %in% brNsign$pos,]
	in_nsign<-tmpNsign[,brdummy$spec[brdummy$group=="in"]]
	if(is.data.frame(in_nsign) && dim(in_nsign)[2]>1){
		mean_inN<-rowMeans(in_nsign,na.rm=T)
	}else if(is.vector(in_nsign) && length(in_nsign)>1){
		mean_inN<-in_nsign
	}else{next}
	if(sum(mean_inN,na.rm=T)==0){next}
	out_nsign<-tmpNsign[,brdummy$spec[brdummy$group=="out"]]
	if(is.data.frame(out_nsign) && dim(out_nsign)[2]>1){
		mean_outN<-rowMeans(out_nsign,na.rm=T)
	}else if(is.vector(out_nsign) && length(out_nsign)>1){
		mean_outN<-out_nsign
	}else{next}
	if(sum(mean_outN,na.rm=T)==0){next}
	ptmp<-wilcox.test(abs(mean_inS-mean_outS),abs(mean_inN-mean_outN))$p.value
	cortex_tab[bpos,2]<-ptmp
	cortex_tab[bpos,3]<-median(abs(mean_inS-mean_outS),na.rm=T) > median(abs(mean_inN-mean_outN),
		na.rm=T)
	signd<-abs(mean_inS-mean_outS)
	nsignd<-abs(mean_inN-mean_outN)
	plot(1,1,type="n",xlim=0:1,ylim=0:1,bty="n",yaxt="n",ylab="",main=b,
		xlab="Openness difference between ingroup and outgroup | Abs (mean ingroup opnenness - mean outgroup openness)")
	axis(2,at=c(.7,.3),labels=c("No acceleration","Acceleration"),tick=F,las=1)
	rect(xleft=quartile_na(signd)[1],xright=quartile_na(signd)[2],ybottom=.2,ytop=.4)
	rect(xleft=quartile_na(nsignd)[1],xright=quartile_na(nsignd)[2],ybottom=.6,ytop=.8)
	segments(y0=c(.3,.7),x0=c(quantile95_na(signd)[1],quantile95_na(nsignd)[1]),
		x1=c(quantile95_na(signd)[2],quantile95_na(nsignd)[2]))
	segments(y0=c(.2,.6),y1=c(.4,.8),x0=c(median_na(signd),median_na(nsignd)),lwd=2)
	text(x=.9,y=.5,signif(ptmp,digits=3))
}
cortex_tab$P<-replace(cortex_tab$P,cortex_tab$P==0,NA)
cortex_tab$Q<-p.adjust(cortex_tab$P,method="BH")
hist(cortex_tab$P,breaks=seq(0,1,.05))
dev.off()

liver<-read.delim("liver_ocrs_filtered_named_hg38_CRE-intersect.bed")
liver$pos<-paste(liver$chr,":",liver$start,"-",liver$end,sep="")
colnames(liver)<-gsub(".","_",colnames(liver),fixed=T)
liverO<-read.delim("liver_ocrs_filtered_named.txt.gz")
colnames(liverO)<-gsub(".","_",colnames(liverO),fixed=T)
liver<-merge(liver,liverO)
liver[6:ncol(liver)]<-replace(liver[6:ncol(liver)],liver[6:ncol(liver)]==-1,NA)

liver_tab<-data.frame(branch=branches,P=double(length(branches)),direction=logical(length(branches)),
	n_sign=double(length(branches)))
pdf("TACIT-liver.pdf",height=4)
for(b in branches){
	bpos<-which(names(branchlist)==b)
	liver_tab[bpos,1]<-b
	brdummy<-data.frame(branch=b,spec=unlist(branchlist[[bpos]]),
		group=c(rep("in",length(unlist(branchlist[[bpos]][1:2]))),
			rep("out",length(unlist(branchlist[[bpos]][3])))))
	brdummy<-brdummy[brdummy$spec %in% colnames(liver),]
	if(nrow(brdummy[brdummy$group=="in",])<1 | nrow(brdummy[brdummy$group=="out",])<1){
		next
	}
	brSign<-sign[sign$branch==b,4:6]
	brNsign<-nsign[nsign$branch==b,c(4,5,9)]
	tmpSign<-liver[liver$pos %in% brSign$pos,]
	in_sign<-tmpSign[,brdummy$spec[brdummy$group=="in"]]
	liver_tab[bpos,4]<-nrow(in_sign)
	if(is.data.frame(in_sign) && dim(in_sign)[2]>1){
		mean_inS<-rowMeans(in_sign,na.rm=T)
	}else if(is.vector(in_sign) && length(in_sign)>1){
		mean_inS<-in_sign
	}else{next}
	if(sum(mean_inS,na.rm=T)==0){next}
	out_sign<-tmpSign[,brdummy$spec[brdummy$group=="out"]]
	if(is.data.frame(out_sign) && dim(out_sign)[2]>1){
		mean_outS<-rowMeans(out_sign,na.rm=T)
	}else if(is.vector(out_sign) && length(out_sign)>1){
		mean_outS<-out_sign
	}else{next}
	if(sum(mean_outS,na.rm=T)==0){next}
	tmpNsign<-liver[liver$pos %in% brNsign$pos,]
	in_nsign<-tmpNsign[,brdummy$spec[brdummy$group=="in"]]
	if(is.data.frame(in_nsign) && dim(in_nsign)[2]>1){
		mean_inN<-rowMeans(in_nsign,na.rm=T)
	}else if(is.vector(in_nsign) && length(in_nsign)>1){
		mean_inN<-in_nsign
	}else{next}
	if(sum(mean_inN,na.rm=T)==0){next}
	out_nsign<-tmpNsign[,brdummy$spec[brdummy$group=="out"]]
	if(is.data.frame(out_nsign) && dim(out_nsign)[2]>1){
		mean_outN<-rowMeans(out_nsign,na.rm=T)
	}else if(is.vector(out_nsign) && length(out_nsign)>1){
		mean_outN<-out_nsign
	}else{next}
	if(sum(mean_outN,na.rm=T)==0){next}
	ptmp<-wilcox.test(abs(mean_inS-mean_outS),abs(mean_inN-mean_outN))$p.value
	liver_tab[bpos,2]<-ptmp
	liver_tab[bpos,3]<-median(abs(mean_inS-mean_outS),na.rm=T) > median(abs(mean_inN-mean_outN),
		na.rm=T)
	signd<-abs(mean_inS-mean_outS)
	nsignd<-abs(mean_inN-mean_outN)
	plot(1,1,type="n",xlim=0:1,ylim=0:1,bty="n",yaxt="n",ylab="",main=b,
		xlab="Openness difference between ingroup and outgroup | Abs (mean ingroup opnenness - mean outgroup openness)")
	axis(2,at=c(.7,.3),labels=c("No acceleration","Acceleration"),tick=F,las=1)
	rect(xleft=quartile_na(signd)[1],xright=quartile_na(signd)[2],ybottom=.2,ytop=.4)
	rect(xleft=quartile_na(nsignd)[1],xright=quartile_na(nsignd)[2],ybottom=.6,ytop=.8)
	segments(y0=c(.3,.7),x0=c(quantile95_na(signd)[1],quantile95_na(nsignd)[1]),
		x1=c(quantile95_na(signd)[2],quantile95_na(nsignd)[2]))
	segments(y0=c(.2,.6),y1=c(.4,.8),x0=c(median_na(signd),median_na(nsignd)),lwd=2)
	text(x=.9,y=.5,signif(ptmp,digits=3))
}
liver_tab$P<-replace(liver_tab$P,liver_tab$P==0,NA)
liver_tab$Q<-p.adjust(liver_tab$P,method="BH")
hist(liver_tab$P,breaks=seq(0,1,.05))
dev.off()

pv<-read.delim("pv_ocrs_filtered_named_hg38_CRE-intersect.bed")
pv$pos<-paste(pv$chr,":",pv$start,"-",pv$end,sep="")
colnames(pv)<-gsub(".","_",colnames(pv),fixed=T)
pvO<-read.delim("pv_ocrs_filtered_named.txt.gz")
colnames(pvO)<-gsub(".","_",colnames(pvO),fixed=T)
pv<-merge(pv,pvO)
pv[6:ncol(pv)]<-replace(pv[6:ncol(pv)],pv[6:ncol(pv)]==-1,NA)

pv_tab<-data.frame(branch=branches,P=double(length(branches)),direction=logical(length(branches)),
	n_sign=double(length(branches)))
pdf("TACIT-PV.pdf",height=4)
for(b in branches){
	bpos<-which(names(branchlist)==b)
	pv_tab[bpos,1]<-b
	brdummy<-data.frame(branch=b,spec=unlist(branchlist[[bpos]]),
		group=c(rep("in",length(unlist(branchlist[[bpos]][1:2]))),
			rep("out",length(unlist(branchlist[[bpos]][3])))))
	brdummy<-brdummy[brdummy$spec %in% colnames(pv),]
	if(nrow(brdummy[brdummy$group=="in",])<1 | nrow(brdummy[brdummy$group=="out",])<1){
		next
	}
	brSign<-sign[sign$branch==b,4:6]
	brNsign<-nsign[nsign$branch==b,c(4,5,9)]
	tmpSign<-pv[pv$pos %in% brSign$pos,]
	in_sign<-tmpSign[,brdummy$spec[brdummy$group=="in"]]
	pv_tab[bpos,4]<-nrow(in_sign)
	if(is.data.frame(in_sign) && dim(in_sign)[2]>1){
		mean_inS<-rowMeans(in_sign,na.rm=T)
	}else if(is.vector(in_sign) && length(in_sign)>1){
		mean_inS<-in_sign
	}else{next}
	if(sum(mean_inS,na.rm=T)==0){next}
	out_sign<-tmpSign[,brdummy$spec[brdummy$group=="out"]]
	if(is.data.frame(out_sign) && dim(out_sign)[2]>1){
		mean_outS<-rowMeans(out_sign,na.rm=T)
	}else if(is.vector(out_sign) && length(out_sign)>1){
		mean_outS<-out_sign
	}else{next}
	if(sum(mean_outS,na.rm=T)==0){next}
	tmpNsign<-pv[pv$pos %in% brNsign$pos,]
	in_nsign<-tmpNsign[,brdummy$spec[brdummy$group=="in"]]
	if(is.data.frame(in_nsign) && dim(in_nsign)[2]>1){
		mean_inN<-rowMeans(in_nsign,na.rm=T)
	}else if(is.vector(in_nsign) && length(in_nsign)>1){
		mean_inN<-in_nsign
	}else{next}
	if(sum(mean_inN,na.rm=T)==0){next}
	out_nsign<-tmpNsign[,brdummy$spec[brdummy$group=="out"]]
	if(is.data.frame(out_nsign) && dim(out_nsign)[2]>1){
		mean_outN<-rowMeans(out_nsign,na.rm=T)
	}else if(is.vector(out_nsign) && length(out_nsign)>1){
		mean_outN<-out_nsign
	}else{next}
	if(sum(mean_outN,na.rm=T)==0){next}
	ptmp<-wilcox.test(abs(mean_inS-mean_outS),abs(mean_inN-mean_outN))$p.value
	pv_tab[bpos,2]<-ptmp
	pv_tab[bpos,3]<-median(abs(mean_inS-mean_outS),na.rm=T) > median(abs(mean_inN-mean_outN),
		na.rm=T)
	signd<-abs(mean_inS-mean_outS)
	nsignd<-abs(mean_inN-mean_outN)
	plot(1,1,type="n",xlim=0:1,ylim=0:1,bty="n",yaxt="n",ylab="",main=b,
		xlab="Openness difference between ingroup and outgroup | Abs (mean ingroup opnenness - mean outgroup openness)")
	axis(2,at=c(.7,.3),labels=c("No acceleration","Acceleration"),tick=F,las=1)
	rect(xleft=quartile_na(signd)[1],xright=quartile_na(signd)[2],ybottom=.2,ytop=.4)
	rect(xleft=quartile_na(nsignd)[1],xright=quartile_na(nsignd)[2],ybottom=.6,ytop=.8)
	segments(y0=c(.3,.7),x0=c(quantile95_na(signd)[1],quantile95_na(nsignd)[1]),
		x1=c(quantile95_na(signd)[2],quantile95_na(nsignd)[2]))
	segments(y0=c(.2,.6),y1=c(.4,.8),x0=c(median_na(signd),median_na(nsignd)),lwd=2)
	text(x=.9,y=.5,signif(ptmp,digits=3))
}
pv_tab$P<-replace(pv_tab$P,pv_tab$P==0,NA)
pv_tab$Q<-p.adjust(pv_tab$P,method="BH")
hist(pv_tab$P,breaks=seq(0,1,.05))
dev.off()

sumTab<-data.frame(sample=c("cortex","liver","pv"),
	isNA=double(3),NS=double(3),agreeing=double(3),disagreeing=double(3),total=rep(84,3))
sumTab[1,2]<-nrow(cortex_tab[is.na(cortex_tab$Q),])
sumTab[1,3]<-nrow(cortex_tab[!is.na(cortex_tab$Q) & cortex_tab$Q>.05,])
sumTab[1,4]<-nrow(cortex_tab[!is.na(cortex_tab$Q) & cortex_tab$Q<.05 & cortex_tab$direction==T,])
sumTab[1,5]<-nrow(cortex_tab[!is.na(cortex_tab$Q) & cortex_tab$Q<.05 & cortex_tab$direction==F,])
sumTab[2,2]<-nrow(liver_tab[is.na(liver_tab$Q),])
sumTab[2,3]<-nrow(liver_tab[!is.na(liver_tab$Q) & liver_tab$Q>.05,])
sumTab[2,4]<-nrow(liver_tab[!is.na(liver_tab$Q) & liver_tab$Q<.05 & liver_tab$direction==T,])
sumTab[2,5]<-nrow(liver_tab[!is.na(liver_tab$Q) & liver_tab$Q<.05 & liver_tab$direction==F,])
sumTab[3,2]<-nrow(pv_tab[is.na(pv_tab$Q),])
sumTab[3,3]<-nrow(pv_tab[!is.na(pv_tab$Q) & pv_tab$Q>.05,])
sumTab[3,4]<-nrow(pv_tab[!is.na(pv_tab$Q) & pv_tab$Q<.05 & pv_tab$direction==T,])
sumTab[3,5]<-nrow(pv_tab[!is.na(pv_tab$Q) & pv_tab$Q<.05 & pv_tab$direction==F,])
write.table(sumTab,"TACIT-summary-table.tsv",quote=F,sep="\t",row.names=F)

plotTab<-t(sumTab[,3:5])
colnames(plotTab)<-c("MotorCortex","Liver","PVInterneuron")
pdf("../../TACIT-barplot.pdf") # Fig. 3C
barplot(prop.table(plotTab,2),horiz=T)
dev.off()

q(save="n")
