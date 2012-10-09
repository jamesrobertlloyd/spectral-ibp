# Load KASP

source('KASP-JL.R')

# Create some synthetic data

n <- 10000
m <- 10
k <- 2
noise <- 0.1

Z <- cbind(matrix(0, n, 1), rbind(matrix(0, n/2, 1), matrix(1, n/2, 1)), rbind(matrix(0, n/4, 1), matrix(1, n/4, 1), matrix(0, n/4, 1), matrix(1, n/4, 1)))
A <- matrix(rnorm((k+1)*m), k+1, m)

X <- Z %*% A
X <- X + noise * rnorm(prod(dim(X)))

# Now run iterative algorithm to learn Z and A from the data

# Initialise variables

Z.hat <- cbind(matrix(1, n, 1), matrix(0, n, k))
A.hat <- matrix(0, k+1, m)

sigma <- 1
alpha <- n / 10

# Estimate A
# determine posterior mean of A instead? sig.x = 1; sig.a = 1; 
# first estimate will be the means: explicitly initialize with means?
lm.hat <- lm(X ~ Z.hat[,2:dim(Z)[2]])
A.hat <- lm.hat$coefficients
A.hat[is.na(A.hat)] <- 0

X.hat <- Z.hat %*% A.hat
X.resid <- X - X.hat

# Now estimate Z1

sp <- KASP(X.resid, sigma=sigma, alpha=alpha)
Z.hat[,2] <- sp - 1

# Re-estimate A

lm.hat <- lm(X ~ Z.hat[,2:dim(Z)[2]])
A.hat <- lm.hat$coefficients
A.hat[is.na(A.hat)] <- 0

X.hat <- Z.hat %*% A.hat
X.resid <- X - X.hat

# Now estimate Z2

sp <- KASP(X.resid, sigma=sigma, alpha=alpha)
Z.hat[,3] <- sp - 1

image(t(Z.hat[seq(1, dim(Z.hat)[1], 100),]))

# Re-estimate A

lm.hat <- lm(X ~ Z.hat[,2:dim(Z)[2]])
A.hat <- lm.hat$coefficients
A.hat[is.na(A.hat)] <- 0

X.hat <- Z.hat %*% A.hat
X.resid <- X - X.hat

# Compute variance explained

print(100 * 1 - (sd(as.vector(X.resid)) / sd(as.vector(X))))

