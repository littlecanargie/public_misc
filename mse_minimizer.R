## Code to verify the minimizer the weighted MSE function does not depend on the weight matrix

# Define loss function (MSE with arbitrary positive definite weight matrix A)
loss <- function(w, A, sig1, sig2, b1, b2){
  W <- matrix(w, 2, 2)
  sum(diag(A %*% (W %*% sig1 %*% t(W) + (diag(2) - W) %*% sig2 %*% t(diag(2) - W) ))) + 
    t(b2 + W %*% (b1 - b2)) %*% A %*% (b2 + W %*% (b1 - b2))
}

# Define function that generates random symmetric positive definite matrices
genposdef <- function(s){
  temp <- matrix(rnorm(s*s),s,s)
  temp %*% t(temp)
}

A <- genposdef(2)
sig1x <- genposdef(2)
sig2x <- genposdef(2)
b1 <- rnorm(2)
b2 <- rnorm(2)

# Check if A is ill-conditioned (eigenvalue close to 0)
eigen(A)$values

# Directly minimizing the loss
round(matrix(optim(c(1, 1, 1, 1), loss, A=A, sig1 = sig1, sig2 = sig2, b1 = b1, b2 = b2)$par, 2, 2), 4)

# The derived formula
round((sig2 + b2 %*% t(b2 - b1)) %*% solve(sig1 + sig2 + (b1 - b2) %*% t(b1 - b2)), 4)