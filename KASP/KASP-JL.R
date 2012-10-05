require(MASS)
require(gregmisc)

#############################################################################
#To implement the Normalized Cut algorithm
#The NCut is equivalent to solving the following generalized eigen system
#	(D-W)y=rDy
#with the constraints that y(i) \in {1,-b} and y'D1=0, where the constraint
#y'D1=0 is automatically satisfied for the generalized eigen system.
#We shall follow the notations used in Shi and Malik
#Mar 26th, 2008
#############################################################################
ncut<-function(W,ncluster,n)
{

#Number of remaining points to NCut
nr<-n;
#Vector to maintain the indices of NCutted points
rvec<-(1:n);
#The similarity matrix left after each round of NCut
rsmat<-matrix(0,n,n);
rsmat<-W;
cut<-NULL;

iloop<-0;
sp<-NA;
while(iloop<ncluster)
{
	iloop<-iloop+1;

	#For the last piece, no further NCut will be performed
	if(iloop==ncluster) 
		{ #cat("iloop =",iloop, "\nThis cut = \n",rvec,"\n\n");
		  sp[rvec]<-iloop;
			next;}

	#Prepare for generalized eigen analysis
	tmp<-rep(1,nr); 
	if(iloop!=1) { rsmat<-rsmat[-cut,]; rsmat<-rsmat[,-cut]; 
	}
	tmp<-rsmat %*% tmp;
	D<-matrix(0,nr,nr); 
	if(length(tmp)!=nr){
		cat("The length of tmp is ",length(tmp),"\n");
		cat("The length of diag(D) is ",nr,"\n");
		browser();
		cat("The length of last cut is ",length(cut),"\n");
		cat("The #rows of rsmat is ",length(rsmat[,1]),"\n");
		cat("The #cols of rsmat is ",length(rsmat[1,]),"\n");
	}

	A<--rsmat; diag(A)<-tmp+diag(A);
	tmp<-1/(sqrt(tmp));
	D<-kronecker(tmp,t(tmp),FUN="*");
	B<-matrix(0,nr,nr);diag(B)<-rep(1,nr);
	A<-A * D;


	#Eigen analysis
	eigens<-eigen(A,symmetric=T);
	ssvector<-eigens$vectors[,nr-1];
	cut<-(1:nr)[ssvector>=0];

	if(F){
	browser();
	z<-matrix(0,1,n);
	tmpp<-1-eigens$values;
	tmpp<-sort(tmpp,decreasing=TRUE);
	for(ii in 1:n-2)
	{ z[ii]<-ii*sum(tmpp[(ii+1):(n-1)]);
	} 
	plot(z[1,1:(n-1)]);
	}


	#We favor the smaller cut
	if(length(cut)>0.5*nr) cut<-(1:nr)[-cut];
	ccut<-rvec[cut];

	#Bad things happen, discard the current loop
	if((length(cut)==n) || (length(cut)==0))
	{ sp[ccut]<-iloop;discard<-1;cat("Discard the current loop\n");break;}

	#To trace back to the original indices of all points currently NCutted
	rvec<-rvec[-cut];
	nr<-nr-length(cut);
	sp[ccut]<-iloop;
	if(nr <1) {cat("Defective sampling\n");break;}
}#End of while(iloop)
sp
}

###########################################################################
# To compute the clustering accuracy using the labels come with the 
# original data. 
# 
# If # clusters <=7, look for maximum match w.r.t. all permutations of 
# 	cluster IDs
# Else maximum match over 10000 random sampling from all permutations
###########################################################################
cRate<-function(sp0, sp1, nc, N)
{
tr<-0;
seqs=seq(1,nc);
spx=matrix(0,N,1);
spy=matrix(0,N,1);

if(nc <8) 
{ perms<-permutations(n=nc,r=nc);
  	np=dim(perms)[1];
  	for(i in 1:np)
  	{
  		for(j in 1:nc) { spx[sp1==j]<-perms[i,j]; }
		tmp<-sum(sp0==spx)/N; 
		if(tr<tmp) {tr=tmp;spy=spx;}
	}
} else
{
  	for(i in 1:10000)
  	{ permx=sample(seqs,nc, replace=FALSE);
  		for(j in 1:nc) { spx[sp1==j]<-permx[j]; }
		tmp<-sum(as.integer(sp0)==spx)/N; if(tr<tmp) {tr=tmp;spy=spx;}
	}
}
cat("The rate is ",tr,"\n");
}

###########################################################################
# alpha is the reduction factor s.t. #(reduced set)=[N/alpha]
KASP <- function(x,ncluster=2,alpha=6000,sigma=0.6)
{
  require(MASS)
  require(gregmisc)
  N<-nrow(x);	#The # of observations
  m<-ncol(x);		#The # of features
  infty<-10^12; # Not sure what thyis is for

  #sp0<-sp;
  y<-x;
  ##############################################################################
  # Start for e-Spectral
  #
  # Implementation based on grouping by Kmeans
  # where cluster centers from Kmeans form a reduced set for fast computing
  # To be extended to k-nearest neighbor ball
  ##############################################################################
  cat("*  # Data points = ",N,"\n");
  cat("*  # Features =",m,"\n");

  n<-floor(length(y[,1])/alpha);
  cat("*  # Representative data points =",n,", data reduction ratio=1/",alpha,"\n");
  cat("*********************************************************************\n");
  cat("Finished data loading......\n");
  cat("Start of KASP @ ",date(),"\n");

  ##############################################################################
  # Run two-stage K-means if there are more than 20000 points
  # else run single K-means
  ##############################################################################
  if(N>20000)
  {
  ############################################################
  #The sampling ratio for the 1st round of K-means
  #0.05 for dataset size larger than 300000
  #this ratio can be even smaller for still larger dataset
  ############################################################
  if(N<300000) { n1<-N*0.1;} else { n1<-N*0.05;}
  idx<-sample((1:N),n1,replace=FALSE);
  xx<-y[idx,];
  cat("Start of coarse K-means ",date(),"\n");
  ###########################################################################
  #The maximum number of iterations and the number of re-starts can both be 
  #adjusted depending on the data. Same applies to other invocations of 
  #K-means 
  ###########################################################################
  xxkms<-kmeans(xx[,1:m],centers=n, iter.max = 200, nstart = 20,
              algorithm = c("Hartigan-Wong"));

  cat("Start of K-means @ ",date(),"\n");
  xkms<-kmeans(y[,1:m],centers=xxkms$centers, iter.max = 200, nstart = 1,
              algorithm = c("Hartigan-Wong"));
  }
  else
  {
  cat("Start of K-means @ ",date(),"\n");
  xkms<-kmeans(y[,1:m],centers=n, iter.max = 200, nstart = 20, algorithm = c("Hartigan-Wong"));
              
  }

  tmp<-xkms$cluster;

  x=xkms$centers;
  smat<-matrix(0,n,n);
  for(i in 1:m)
  {
          tmp<-kronecker(x[,i],t(x[,i]),FUN="-");
          tmp1<-matrix(tmp,n,n);
          smat<-smat+tmp1*tmp1;
  }
  #}
  smatz<-smat;

  ######################################################################
  # The bandwidth parameter (sigma)
  # This value is to be searched over a range
  # The following, in the form of (reduction ratio, sigma), are examples 
  # of good (by no means the best) values based on our experiments.
  # Note even for the same dataset, the optimal sigma varies for different 
  # reduction ratios
  #
  # Musk			(1/8,30)
  # Magic Gamma		(1/8,36)
  # Connect-4		(1/200,50) 
  # USCI			(1/500,10)
  # Poker Hand		(1/3000,0.6)
  ######################################################################
  #sigma<-0.6;
  cat("Sigma=",sigma,"\n");
  smat<-matrix(0,n,n);
  smat<-exp(-smatz/sigma);
  tmp<-NULL;tmp1<-NULL;

  sp<-ncut(smat,ncluster,n);
  #To recover the cluster membership for all data points
  spp<-xkms$cluster;
  for(i in 1:n)
  {
	  spp[spp==i]<-sp[i];
  }
  cat("End of KASP @ ",date(),"\n");

  return(spp) # Return the cluster assignments
  #cRate(sp0,spp,ncluster,N);
  ####}#End of Loop on sigma
}
