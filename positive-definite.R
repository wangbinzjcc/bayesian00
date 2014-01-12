
##########################################################
# Generating a random positive-definite matrix with user-specified positive
# eigenvalues
# If eigenvalues are not specified, they are generated from a uniform
# distribution
########################################################### 
Posdef <- function(n, ev = runif(n, 0, 10)){ 
  Z <- matrix(ncol=n, rnorm(n^2)) 
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp) 
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph) 
  Z <- t(O) %*% diag(ev) %*% O
  return(Z) } 
pdmat <- Posdef(n=3, ev=3:1) 
eigen(pdmat)$val
###


