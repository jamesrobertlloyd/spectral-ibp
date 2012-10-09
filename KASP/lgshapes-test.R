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
  plot.A(A)
  # select the latent features for each datum
  Z = matrix(runif(N*K) > Zprob, nrow=N)
  
  # generate the data
  gen.norm <- function(mean, sd, n=1){ rnorm(n, mean, sd) }
  means = Z %*% A
  X = apply(means, 1:2, FUN=gen.norm, sd=sigX)
  
  return(list(X=X, A=A, Z=Z, K=K, N=N, D=D, sigX=sigX))
}

# load needed packages
source('KASP-JL.R')

# initialization
model = gen.data(1000, sigX=1)
X = model$X
params = list(K.init=4, sigX=1.0, sigA=1.0, n.iter=100, kasp.sig=1.0, kasp.alpha=model$N/100)
Z.hat = cbind(rep(1, N), matrix(0, model$N, params$K.init - 1))
A.hat = map.A(X, Z.hat, params$sigX, params$sigA) 
X.hat = Z.hat %*% A.hat
X.resid = X - X.hat
Kval = params$K.init

# latent feature assignment iterations
x11()
for(itn in 1:params$n.iter){
  for (k in 1:Kval){
    if(itn == 1 && k == 1){ next }
   
    sp <- KASP(X.resid, sigma=params$kasp.sig, alpha=params$kasp.alpha)
    Z.hat[,k] <- sp - 1
    A.hat = map.A(X, Z.hat, params$sigX, params$sigA)
    plot.A(A.hat)
    rm.ind = ifelse(k == Kval, 1, k + 1)
    Z.hat[,rm.ind] = 0
    X.hat = Z.hat %*% A.hat
    X.resid = X - X.hat
  }
}




