map.A <- function(X, Z.hat, sigX, sigA){
  solve(t(Z.hat) %*% Z.hat + (sigX/sigA)^2*diag(dim(Z.hat)[2]))%*%t(Z.hat)%*%X
}

plot.A <- function(mat){
  nfeat = dim(mat)[1]
  rdim = sqrt(dim(mat)[2])
  nr = ceiling(sqrt(nfeat))
  nc = floor(sqrt(nfeat))
  par(mfrow=c(nr, nc))
  for (i in 1:nfeat){
    image(t(matrix(mat[i,], nrow=rdim)), axes=FALSE)  
  }
  
}

gen.data <- function(N, K=4, D=36, Zprob = 0.5, sigX=0.1){
  # define four latent features and generate data
  A = matrix(0, K, D)
  feats.inds = list(c(2, 7, 8, 9, 14), c(4, 5, 6, 10, 12, 16, 17, 18), c(19, 25, 26, 31, 32, 33), c(22, 23, 24, 29, 35))
  par(mfrow=c(2,2))
  for (i in 1:length(feats.inds)){
    A[i, feats.inds[[i]]] = 1 
  } 
  # TODO make an option
  # x11()
  # plot.A(A)
  # select the latent features for each datum
  Z = matrix(runif(N*K) > Zprob, nrow=N)
  
  # generate the data
  gen.norm <- function(mean, sd, n=1){ rnorm(n, mean, sd) }
  means = Z %*% A
  X = apply(means, 1:2, FUN=gen.norm, sd=sigX)
  
  return(list(X=X, A=A, Z=Z, K=K, N=N, D=D, sigX=sigX))
}

calc.llhood <- function(X, Z, K, sigX, sigA){
  ZtZ = t(Z) %*% Z
  N = dim(Z)[1]
  D = dim(X)[2]
  diag.sig = diag((sigX/sigA)^2, K)
  -N*D/2*log(2*pi) - (N-K)*D*log(sigX) - K*D*log(sigA) - D/2*log(det(ZtZ+ diag.sig)) - 1/(2*sigX^2)*sum(diag(t(X) %*% (diag(N) - Z %*% solve(ZtZ + diag.sig) %*% t(Z)) %*% X))
}

# load needed packages
source('KASP-JL.R')

# initialization
trainF = 4/5
model = gen.data(N=1000*1/trainF, sigX=1.0)
trN = ceiling(model$N*4/5)
X = model$X[1:trN,]
Xtest = model$X[(trN + 1):model$N,]
params = list(K.init=4, sigX=1.0, sigA=1.0, n.iter=20, kasp.sig=1.0, kasp.alpha=trN/11)
Z.hat = cbind(rep(1, trN), matrix(0, trN, params$K.init - 1))
A.hat = map.A(X, Z.hat, params$sigX, params$sigA) 
X.hat = Z.hat %*% A.hat
X.resid = X - X.hat
Kval = params$K.init

# init for latent feature assignments
llhood = c()
max.llhood = -Inf
x11()
dev.A = as.integer(dev.cur())
x11()
par(mfrow=c(1,1))
dev.llhood = as.integer(dev.cur())

# latent feature assignment iterations
for(itn in 1:params$n.iter){
  for (k in 1:Kval){
    if(itn == 1 && k == 1){ next }
    # perform spectral clustering on residual matrix
    sp <- KASP(X.resid, sigma=params$kasp.sig, alpha=params$kasp.alpha)
    Z.hat[,k] <- sp - 1
    A.hat = map.A(X, Z.hat, params$sigX, params$sigA)
    
    # display updated latent features
    dev.set(dev.A)
    plot.A(A.hat)
    
    # likelihood calculations/display
    llh = calc.llhood(X, Z.hat, Kval, params$sigX, params$sigA)
    llhood[length(llhood) + 1] = llh
    dev.set(dev.llhood)
    plot(llhood,type='b',lwd=2)
    if (llh > max.llhood){
      Z.hatML = Z.hat
      max.llhood = llh
    }
    # remove the next Z column and calculate the residuals
    rm.ind = ifelse(k == Kval, 1, k + 1)
    Z.hat[,rm.ind] = 0
    X.hat = Z.hat %*% A.hat
    X.resid = X - X.hat
  }

}
A.hatML = map.A(X, Z.hatML, params$sigX, params$sigA)
plot.A(A.hatML)



