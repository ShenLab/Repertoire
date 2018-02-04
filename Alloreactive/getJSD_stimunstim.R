rm(list=ls())
options(stringsAsFactors=FALSE)

# sample name and subtype of interest
L=args[1]
type='CD4'
cutoff=NA # top N clones. Set to NA if no cutoff desired

# load file and process data
D=read.table(L,header=T)
if(is.na(cutoff)){
    DD=data.frame(unstim=D[[paste(type,'fq.unstim',sep='')]],stim=D[[paste(type,'fq.stim',sep='')]])
    rownames(DD)=D$nucleotide
    DD=subset(DD, unstim>0 | stim>0) # remove 0 entries corresponding to other subtype information
} else {
    # extract stim and unstim separately
    DDunstim=data.frame(nucleotide=D$nucleotide,unstim=D[[paste(type,'fq.unstim',sep='')]])
    DDstim=data.frame(nucleotide=D$nucleotide,stim=D[[paste(type,'fq.stim',sep='')]])       
    
    # get top N clones as defined by cutoff
    un=DDunstim[order(DDunstim$unstim,decreasing=T),]
    st=DDstim[order(DDstim$stim,decreasing=T),]
    un=un[1:max(which(un$unstim==un$unstim[cutoff])),]
    st=st[1:max(which(st$stim==st$stim[cutoff])),]
    
    # merge stim and unstim
    DD=merge(un,st,by='nucleotide',all=T) # merge stim and unstim
    DD[is.na(DD)]=0
    DD=subset(DD, unstim>0 | stim>0) # remove 0 entries corresponding to other subtype information
    rownames(DD)=DD$nucleotide
    DD=DD[,c(2,3)]
}

DD=sweep(DD, 2, colSums(DD), FUN="/")  # make sure everything is normalized 

# define entropy
entropy<-function(p)
{
  p=p/sum(p) # normalization, probably redundant here
  p=p[p>0] # return 0s
  return(sum(-p*log2(p)))
}

# combine into a new distribution
M=(DD$unstim+DD$stim)/2

# Jenson Shannon Divergence (add a square root for distance)
JSD=entropy(M)-0.5*(entropy(DD$unstim)+entropy(DD$stim))
print(JSD)