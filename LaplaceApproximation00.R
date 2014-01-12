##### ====================================================================
require(LaplacesDemon)
### ==================================

### =====================================
Model <- function(parm, Data)
{  beta <- parm[1:Data$J]
  sigma <- exp(parm[Data$J+1])
  ### 
  beta.prior <- dnormv(beta, 0, 1000, log=TRUE)
  sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
  ### 
  mu <- tcrossprod(beta, Data$X)
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### 
  LP <- LL + sum(beta.prior) + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}
### =======================================
N <- 100
J <- 5
X <- matrix(1, N, J)
for (j in 2:J) {X[, j] <- rnorm(N, runif(1, -3, 3), runif(1, 0.1, 1))}
beta <- runif(J, -3, 3)
e <- rnorm(N, 0, 0.1)
y <- tcrossprod(X, t(beta)) + e
# 
#
#
mon.names <- c("LP", "sigma")
parm.names <- as.parm.names(list(beta=rep(0, J), log.sigma=0))
PGF <- function(Data) return(c(rnormv(Data$J, 0, 10), log(rhalfcauchy(1, 25))))
### 
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
               parm.names=parm.names, y=y)
### =========================================================
### =========================================================
    Iterations=100
    Stop.Tolerance=1.0E-4
### 
GIV=function (Model, Data, n = 100 ){
  LIV <- length(Data$parm.names)
  iv <- rep(NA, LIV)
  for (i in 1:n) {# i=1 
    IV <- Data$PGF(Data)
    M <- try(Model(IV, Data), silent = TRUE)
    if (!inherits(M, "try-error") & is.finite(M[[1]]) & is.finite(M[[2]]) & 
          identical(as.vector(M[[5]]), IV)) 
    {   iv <- IV
        break
    }          }
  if ((i == n) | any(is.na(iv))) 
    cat("\nWARNING: Acceptable initial values were not generated.\n")
  return(iv)                        }

#

Initial.Values <- GIV(Model, Data, PGF=TRUE)
     m.old <- Model(Initial.Values, Data)
     parm <- m.old[["parm"]]
 #
######### Begin Laplace Approximation #############
 LBFGS <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     { 
       Dev <- matrix(m.old[["Dev"]],1,1)
       parm.len <- length(as.vector(parm))
       parm.new <- parm.old <- m.old[["parm"]]
       names(parm.new) <- names(parm.old) <- Data$parm.names
       tol <- 1
       post <- matrix(parm.new, 1, parm.len)
  ModelWrapper <- function(parm.new) {
       out <- Model(parm.new, Data)[["LP"]]
       return(out)                   }
   #    
    for (iter in 1:Iterations) {
        parm.old <- parm.new
        ### LBFGS
        Fit <- optim(par=parm.new, fn=ModelWrapper,
                      method="L-BFGS-B", control=list(fnscale=-1, maxit=1)
                      )
        m.new <- Model(Fit$par, Data)
    if(!is.finite(m.new[["LP"]]) | {m.new[["LP"]] < m.old[["LP"]]})
       { m.new <- m.old }
       parm.new <- m.new[["parm"]]
       post <- rbind(post, parm.new)
       Dev <- rbind(Dev, m.new[["Dev"]])
       tol <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
      if(tol <= Stop.Tolerance) break
                              }
      Dev <- Dev[-1,]; post <- post[-1,]
       ### Output
      LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
                  parm.old=parm.old, post=post, Step.Size=0,
                  tol.new=tol)
       return(LA) 
     }
### ==========================================================