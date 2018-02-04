#########
#
# Repertoire fitting for T cell repertoire to convert read counts into template counts
#
# To run: Rscript RepertoireFitting.R <full path to sample file> <output path for plotting figures> <type> <verbosity>
#
#  INPUT
#   type: column name for data
#   verbosity: 0 -> shows only final result, 1 -> result at each step of optimization
#
#  OUTPUT
#    minimum KL divergence
#    parameters: depth, overdispersion, power law slope
#
#
#   EXAMPLE: Rscript reads2templates.R myfile tmp/ CD4.unstim 1
#
#########

###### Functions to generate overdispersed reads given parameters coverage, degree of overdispersion, and powerlaw slope 

rm(list=ls())
options(stringsAsFactors=FALSE)

args <- commandArgs(TRUE)
f=args[1] # path to file
figurepath <-args[2] #  args[2]  # plotting path 
type=args[3]
verbosity=as.integer(args[4])
nruns=5
iter=5000 # max iterations

###### Begin simulation functions ######
# compute the initial powerlope slope (a)
geta <- function(distr,depth)
{
  # Estimate number of cells/templates
	distr_templates=distr
#	distr_templates=distr_templates[distr_templates$counts>2,] # remove counts of 2
	counts=ceiling(distr_templates$counts/depth)
	distr_templates$counts=counts
	distr_templates=aggregate(Y~counts,distr_templates,sum)

	# select the power law subset by frequency
	ind=min(which(distr_templates$Y==1))
	distr_templates=distr_templates[seq(ind),]
	
	# linear fit
	fit=lm(log(distr_templates$Y)~log(distr_templates$counts))
	slope=coefficients(fit)[[2]]
	
	return(-slope)
}

# simulate overdispersed reads
simulatereads<-function(x,coverage,sp,a)
{
  # get number of clones for each cell count
  p=x^(-a)
  p=p/sum(p)
  count=round(p*nclones)
  count=count[count>0]
  xx=seq(count)
  
  # negative binomial sampling
      # size is linearly dependent on read counts
      # mean is the read counts (cell counts * coverage)
  r=unlist(lapply(xx, function(i) rnbinom(count[i],size=sp,mu=(i*coverage))))

  return(r)
}

# histogram of counts for abundance distr
counts2hist<-function(c)
{
  c=c[c>0]
  counts=c
  c=c
  h<-table(c)
  h<-cbind(names(h),h)
  mode(h) <- "numeric"
  h<-as.data.frame(h)
  names(h)<-c("counts","Y")
  return(h)
}
###### End simulation functions ######

##### Functions to compute information content (Shannon entropy, Kullback-Leibler) #####

# function to compute entropy
calcH<-function(p)
{
	p=p/sum(p)
	H=-sum(p*log(p))
	return(H)
}

# function to compute KL divergence
calcL<-function(p,q)
{
      # normalize
      p=p/sum(p)
      q=q/sum(q)

      Hp=p*log(p)
      Hpq=p*log(q)
      f=Hp-Hpq
      f[is.nan(f)]=0
      
    return(sum(f))
}

# function to preprocess the data in order to compute KL divergence
getL<-function(X,Y=aa)
{  
  aa=Y # real data distr
  bb=counts2hist(X) # distr of simulated data
 
  # Initialize arrays. Each index (1,2,3,...) is a cell copy number, starting at 2 as in real data.
  xmax=max(aa$counts,bb$counts)
  Mp=array(0,dim=(xmax-1))
  Mq=array(0,dim=(xmax-1))
 
  # Populate arrays with clone copy numbers 
  for(i in seq(dim(bb)[1])) {Mq[bb$counts[i]-1]=bb$Y[i]} # simulated
  for(i in seq(dim(aa)[1])) {Mp[aa$counts[i]-1]=aa$Y[i]} # real

 
  # remove q=0 entries to avoid issues in KL computation 
  inds=which(Mq!=0)
  Mp=Mp[inds]
  Mq=Mq[inds]

  # compute the divergence
  L=calcL(Mp,Mq)
 
  return(L)
}
###### End analysis functions ######


###### Optimization #######

# Updating parameters (slope, coverage, overdispersion) via random walk
paramupdate<-function(theta,i,stepsize)
{
    direction=round(runif(1))*2-1 # pick a direction
    step=stepsize[i] # get step size
    if ((theta[i]+step*direction)>10^-10){ # check that parameters are greater than zero. Also deals with numeric issues.
      theta[i]=theta[i]+step*direction
    }
    return(theta)
}

# Simulated Annealing
nsteps=1 # Keep track of steps taken with the same parameter settings (convergence criteria)
runMCMC<-function(x_init,theta_init,slope,iter=1000,temp=100,Lratio=0.75,stepsize)
{
  x=x_init
  # initial parameters
  nparams = length(theta_init)
  chain = array(dim=c(iter+1,nparams+4))  ## CHAIN <- Variable stores information on each run of the optimization
						# 5 columns (iteration, coverage, overdispersion, slope, KL)
  theta=theta_init;
  r=simulatereads(x=x,coverage=theta[1],sp=theta[2],a=slope)
  Lupper=getL(r)
  L=Lupper;
  chain[1,] = c(0,theta_init,slope,L,temp) # Starting parameters at time 0
 
  if(verbosity>=1){
    cat('iter\tdepth\tod\tslope\tKL\ttemp\n')
    cat(chain[1,],'\n')
      }
  
  # iteratively update parameters
  for(i in seq(iter)){
    thetanew=paramupdate(theta,i%%2+1,stepsize)

    if(sum(theta-thetanew)==0){ # do not update if parameters are the same (i.e. bad parameter space)
      next;
    }
    
    # simulate data and compute likelihood
    r=simulatereads(x=x,coverage=thetanew[1],sp=thetanew[2],a=slope)
    slope=geta(aa,theta[1])
    Lnew=getL(r)
    tol=1e-4
    if ((Lnew<L) & abs(L-Lnew)>tol){ # take the step if divergence decreases
    	L=Lnew
  	  theta=thetanew
    	nsteps=1
      }
    else if((L<Lnew) & abs(L-Lnew)>tol ){ # decide whether to take a step
	if(runif(1)<(exp(temp*(L-Lnew)))){
	  L=Lnew
	  theta=thetanew
	  nsteps=1
	} else{nsteps=nsteps+1} # if parameters haven't changed, add 1
    }
    
    # check for convergence after a burnin period
    if(nsteps>100 & i>1000) {return(chain)}
    
    # cooling step
    if((L/Lupper)<Lratio){
      temp=temp*5 # lower probability threshhold
      stepsize=stepsize/2
      Lupper=L}
    
    # add to the chain
    chain[i+1,] = c(i,theta,slope,L,temp)
  	if(verbosity>=1){cat(c(i,theta,slope,L,temp,'\n'))}
    }
      	return(chain)
}

#### End optimization #### 


########### RUN THE MODEL ###########
###### Load real data and initialize data dependent parameters ###### 
R<-read.table(f,header=T)
aa<-counts2hist(R[[type]])
nR=sum(R[[type]])
nclones=sum(aa$Y)

###### End initialization ######

sp_init=1 # overdispersion
covvals=c(10,50,100,200,ceiling(nR/nclones)) # initial guess

bestL=Inf
for(j in seq(nruns)){
  coverage_init=covvals[j]
  a_init=geta(aa,coverage_init)
  theta0=c(coverage_init,sp_init)
  besttheta=rep(0,length(theta0))
  nmat=array(dim=c(nruns,length(theta0)+4))
  
output=runMCMC(seq(5000),theta0,a_init,iter=iter,stepsize=c(coverage_init/10,0.1),temp=5)  # THIS LINE RUNS THE OPTIMIZATION
converged=output[max(which(!is.na(rowSums(output)))),]

nmat[j,]=converged
L=converged[5]
    if(L<bestL){
      bestL=L
      besttheta=converged[c(2,3,4)]
    }
}

cat('minKL = ' , bestL, '\n\n')
cat('depth, od, slope\n')
cat(besttheta[1],',',besttheta[2],',',besttheta[3],'\n\n')


pdf(paste(figurepath,'/',f,'.pdf',sep=''),pointsize=12)
r=simulatereads(x=seq(5000),coverage=besttheta[1],sp=besttheta[2],a=besttheta[3])
bb=counts2hist(ceiling(r/besttheta[1]))
ntemplates=sum(bb$Y*bb$counts)
par(mar=c(5,5,2,2))
plot(bb$counts,bb$Y,log='xy',xlim=c(1,max(c(R[[type]],r))),xlab='copy number', ylab='clones',cex=2,cex.main=2,cex.lab=2,cex.axis=2,pch=20)
out<-dev.off()
