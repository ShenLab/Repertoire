rm(list=ls())
options(stringsAsFactors=FALSE)
args <- commandArgs(TRUE)


type="CD8.stim"
L=args[1]
S=read.table(L,header=T) # load up all samples

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
  slope=coefficients(fit)[[2]]
  return(slope)
}

calcH<-function(vals){
  vals=vals[vals>0]
  fq=vals/sum(vals)
  H=-sum(fq*log2(fq))
  return(H)
}

calcC<-function(vals){
  vals=vals[vals>0]
  H=calcH(vals)
  Hmax=log2(length(vals))
  C=1-H/Hmax
  return(C)
}

calcSI<-function(vals){
  vals=vals[vals>0]
  fq=vals/sum(vals)
  si=sum(fq^2)
  return(si)
}

calcr20 = function(X){
  X=sort(X,decreasing=T)
  X=X[X>0]
  CX=cumsum(X)
  num=length(which(CX/sum(X)<=0.2))
  den=length(X)
  return(num/den)
}

calcr50 = function(X){
  X=sort(X,decreasing=T)
  X=X[X>0]
  CX=cumsum(X)
  num=length(which(CX/sum(X)<=0.5))
  den=length(X)
  return(num/den)
}

clipboard <- function(x, sep="\t", row.names=FALSE, col.names=TRUE){
  con <- pipe("pbcopy -selection clipboard -i", open="w")
  write.table(x, con, sep=sep, row.names=row.names, col.names=col.names)
  close(con)
}

count=S[[paste(type,sep='')]]
count=count[count>0]
name=strsplit(basename(L),split="[.]")[[1]][1]

summarystats=data.frame(samples=name)
summarystats$maxfq=signif(max(count)/sum(count),2) # max frequency
summarystats$entropy=signif(calcH(count),2) # entropy
summarystats$clonality=signif(calcC(count),2) # clonality
summarystats$SI=signif(calcSI(count),2) # simposn index
summarystats$r20=signif(calcr20(count),2) #r20
summarystats$r50=signif(calcr50(count),2) #r50
summarystats$slope=getslope(counts2hist(count)) #slope
clipboard(summarystats,col.names=F) # save to clipboard
print(summarystats)