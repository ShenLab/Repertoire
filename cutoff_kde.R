rm(list=ls())
options(stringsAsFactors=FALSE)
args <- commandArgs(TRUE)

library(MASS)
f='YOURFILEHERE'
name=strsplit(f,"[.]")[[1]][1]
type='CD8'
pseudocount=0.5

mydata=read.table(f,header=T)
mydata=data.frame(nucleotide=mydata$nucleotide, unstim=mydata[[paste(type,'.unstim',sep="")]], stim=mydata[[paste(type,'.stim',sep="")]])
mydata$unstimfq=mydata$unstim/sum(mydata$unstim)
mydata$stimfq=mydata$stim/sum(mydata$stim)
originaldata=mydata

mydata[mydata==0]=pseudocount
mydata$unstimfq=mydata$unstim/sum(mydata$unstim)
mydata$stimfq=mydata$stim/sum(mydata$stim)


xr=range(mydata$unstimfq)
yr=range(mydata$stimfq)

sharedclones=subset(mydata,stim>pseudocount & unstim>pseudocount)
sharedoriginal=subset(originaldata,stim>0 & unstim>0)

uniqueclones=subset(mydata,stim==pseudocount | unstim==pseudocount)

par(mar=c(5,5,3,3))

sharedclones$foldchange=sharedclones$stimfq/sharedclones$unstimfq
sharedoriginal$foldchange=sharedoriginal$stimfq/sharedoriginal$unstimfq

tt=subset(sharedclones,foldchange<=Inf)
ttt=subset(tt,stimfq>5e-5 | unstimfq>5e-5)
nucs=ttt$nucleotide

orig=originaldata[which(originaldata$nucleotide %in% nucs),]


yhist=table(signif(tt$foldchange,2))
h<-cbind(names(yhist),yhist)
mode(h) <- "numeric"
h<-as.data.frame(h)
names(h)<-c('x','y')

par(mar=c(5,5,2,2))
ss=data.frame(un=tt$unstimfq,st=tt$stimfq)

originaldata=subset(originaldata,unstimfq>0 & stimfq>0)
smoothScatter(log10(originaldata$unstimfq),log10(originaldata$stimfq),
                    nbin=100,nrpoints=0,bandwidth=c(0.2,0.2), xlim=c(-5.2,-1),ylim=c(-5.2,-1),colramp = colorRampPalette(c("white", "white")),cex=2,xlab='',ylab='',cex.axis=2)

z <- kde2d(log10(orig$unstimfq), log10(orig$stimfq), n=100)
points(log10(orig$unstimfq), log10(orig$stimfq), col = 'black',pch=20,cex=3,cex.axis=2.5)
contour(z, drawlabels=FALSE, nlevels=15, col="red", add=TRUE,lwd=4)

lines(c(-5,-4,-3,-2,-1),1*c(-5,-4,-3,-2,-1),col='purple',lwd=5,lty=2)
lines(c(-5,-4,-3,-2,-1),log10(2*10^(c(-5,-4,-3,-2,-1))),col='purple',lwd=5,lty=1)
lines(c(-5,-4,-3,-2,-1),log10(3*10^(c(-5,-4,-3,-2,-1))),col='purple',lwd=5,lty=2)
lines(c(-5,-4,-3,-2,-1),log10(4*10^(c(-5,-4,-3,-2,-1))),col='purple',lwd=5,lty=2)
lines(c(-5,-4,-3,-2,-1),log10(5*10^(c(-5,-4,-3,-2,-1))),col='purple',lwd=5,lty=2)