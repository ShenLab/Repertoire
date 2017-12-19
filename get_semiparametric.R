rm(list=ls())
options(stringsAsFactors=FALSE)
args <- commandArgs(TRUE)

counts2hist<-function(c)
{
  c=c[c>0]
  counts=c
  c=c/sum(c)
  h<-table(c)
  h2=table(counts)
  h<-cbind(names(h),h)
  h2<-cbind(names(h2),h2)
  mode(h) <- "numeric"
  mode(h2) <- "numeric"
  h<-as.data.frame(h)
  h2<-as.data.frame(h2)
  names(h)<-c("fq","Y")
  names(h2)<-c("counts","Y")
  h$counts=h2$counts
  return(h)
}
getslope<-function(vals)
{
  vals=vals[vals$counts>0,]
  if(length(which(vals$Y==1))==1){
    vals2=vals
  }else{
    if(abs(diff(log(vals$fq[sort(which(vals$Y==1))[1:2]])))<1.5){
      ind=sort(which(vals$Y==1))[2]
    }
    else{
      ind=sort(which(vals$Y==1))[1]
    }
    vals2=vals[1:ind,]
  }
  fq=vals2$fq
  fit=lm(log(vals2$Y)~log(fq))
  return(fit)
}

# Load files
f=args[1]
name=strsplit(basename(f),split="[.]")[[1]][1]
type='CD4'
D=read.table(f,header=T)
unstim=data.frame(nucleotide=D$nucleotide,unstim=D[[paste(type,'.unstim',sep='')]])
stim=data.frame(nucleotide=D$nucleotide,stim=D[[paste(type,'.stim',sep='')]])

# Find shared clones in stim and unstim
Dmerge=merge(stim,unstim,by='nucleotide',all=T)
Dmerge[is.na(Dmerge)]=0
Dmerge$unstimfq=Dmerge$unstim/sum(Dmerge$unstim)
Dmerge$ind=0+(Dmerge$stim>0 & Dmerge$unstim>0)

DD=subset(Dmerge,ind==1)
nseen=dim(DD)[1]

un=counts2hist(unstim$unstim)
st_in_un=counts2hist(DD$unstim)
st_in_un$fq=st_in_un$counts/sum(unstim$unstim)

# get powerlaw parameters
f.un=getslope(un)
f.st=getslope(st_in_un)
s.un=signif(coefficients(f.un)[[2]],2)
s.st=signif(coefficients(f.st)[[2]],2)


# fit powerlaw to entire dataset
newfq=c(1e-7,un$fq)
yy=exp(predict(f.un,newdata=data.frame(fq=newfq)))

newfq_stim=c(1e-7,st_in_un$fq)
yyst=exp(predict(f.st,newdata=data.frame(fq=newfq_stim)))

mydata=merge(stim,unstim,all=T)
data_stim=subset(mydata,stim>0)
observed=subset(data_stim,unstim>0)

# number stim not seen in unstim
nzeros=sum(dim(data_stim)-dim(observed))

# compute the intercept
xintercept=exp((log(nzeros)-coefficients(f.st)[[1]])/coefficients(f.st)[[2]])

###PLOTTING
library(scales)
par(mar=c(5,5,3,2))
plot(un$fq,un$Y,log='xy',pch=16,cex=3,xaxt='n',yaxt='n',xlab='',ylab='',main='',cex.main=2,cex.lab=2,cex.axis=2,xlim=c(5e-8,5e-2),ylim=c(1,1e6))
box(lwd=2)
points(st_in_un$fq,st_in_un$Y,pch=21,col='blue',bg=alpha('lightblue',0.8),cex=3)
lines(newfq,yy,col='red',lwd=8,lty=2)
lines(newfq_stim,yyst,col='red',lwd=8,lty=2)

axis(side = 2, las=2,at=c(1,10,1e2,1e3,1e4,1e5,1e6), labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6)),cex.axis=2)
axis(side = 1, las=1,at=c(1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1), labels = c(expression(10^-7),expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)),cex.axis=2,padj=0.2)

abline(h=nzeros,col='darkgreen',lwd=8,lty=2)
points(xintercept,nzeros,col='darkgreen',pch=16,cex=5,lwd=3)
abline(v=xintercept,lty=3,lwd=8,col='purple')

cat('x-intercept = ',xintercept,',\t','unseen = ',nzeros,'\n')

