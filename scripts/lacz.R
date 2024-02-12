dat<-read.delim("data/lacZ-results-table.tsv")

datb<-aggregate(. ~ Embryo,dat[,-3],max)
datb[,3:10]<-replace(datb[,3:10],datb[,3:10]>1,1)

fisher.p<-double(length=ncol(datb)-2)
for(f in 1:(ncol(datb)-2)){
	fisher.p[f]<-fisher.test(table(datb[,c(2,f+2)]))$p.value
}

notochord<-table(datb[,c(2,3)])
notochord<-t(notochord/rowSums(notochord))
notochord<-notochord[2:1,]
drg<-table(datb[,c(2,4)])
drg<-t(drg/rowSums(drg))
drg<-drg[2:1,]
somites<-table(datb[,c(2,5)])
somites<-t(somites/rowSums(somites))
somites<-somites[2:1,]
fb<-table(datb[,c(2,6)])
fb<-t(fb/rowSums(fb))
fb<-fb[2:1,]
fAER<-table(datb[,c(2,7)])
fAER<-t(fAER/rowSums(fAER))
fAER<-fAER[2:1,]
hAER<-table(datb[,c(2,8)])
hAER<-t(hAER/rowSums(hAER))
hAER<-hAER[2:1,]
CMB<-table(datb[,c(2,9)])
CMB<-t(CMB/rowSums(CMB))
CMB<-CMB[2:1,]
aEye<-table(datb[,c(2,10)])
aEye<-t(aEye/rowSums(aEye))
aEye<-aEye[2:1,]

col<-c(grey(.5),grey(1))
pdf("lacZ-relative.pdf")
par(mfrow=c(2,4),las=1)
barplot(notochord,col=col,main="Neural tube",names.arg=c("ancestral","derived"))
barplot(drg,col=col,main="Dorsal root ganglia",names.arg=c("ancestral","derived"))
barplot(somites,col=col,main="Somites",names.arg=c("ancestral","derived"))
barplot(fb,col=col,main="Forebrain",names.arg=c("ancestral","derived"))
barplot(fAER,col=col,main="Forelimb AER",names.arg=c("ancestral","derived"))
text(1,.9,"P = 0.0027")
barplot(hAER,col=col,main="Hindlimb AER",names.arg=c("ancestral","derived"))
text(1,.9,"P = 0.00018")
barplot(CMB,col=col,main="Common Midbrain Background",names.arg=c("ancestral","derived"))
barplot(aEye,col=col,main="Supra-ocular",names.arg=c("ancestral","derived"))
dev.off()

q(save="no")
