#################################################################
#     weaver 2008 R code            wangbinzjcc 2014-1-13
##################################################################

# APPENDIX INCLUDING WinBUGS AND R CODE
#  WinBUGS Code for Logistic Regression Model Fit to Ingot Experimental Data
# r - number of ingots not ready for rolling
# n - number of tested ingots
# N - number of experimental runs
# pi - probability of ingot not ready for rolling or failure probability
# ST - soaking time
# HT - heating time
# X1 - centered soaking time
# X2 - centered heating time
# 
require(LaplacesDemon)     
#
Model <- function(parm, Data){
# parameters
  beta0 <- parm[1]
 # beta1 <- parm[2]
  beta2 <- parm[2]
  beta12 <- parm[3]
  beta11 <- parm[4]
  beta22 <- parm[5]  
  #
  ST <- Data$ST
  HT <- Data$HT
  r <- Data$r0
  n <- Data$n0
  X1 <- ST - mean(ST)
  X2 <- HT - mean(HT)
# log-likelihood
  logit.pi <- beta0 + #beta1*X1 + 
              beta2*X2 + beta12*X1*X2 + beta11*X1*X1 + beta22*X2*X2
  pi <- invlogit(logit.pi)
  LL <- sum(dbinom(x=r, size=n, prob=pi, log=T))
# priors
  beta0.priors <- dnormv(beta0, 0.0, 1.0E3, log=T)
#  beta1.priors <- dnormv(beta1, 0.0, 1.0E3, log=T)
  beta2.priors <- dnormv(beta2, 0.0, 1.0E3, log=T)
  beta12.priors <- dnormv(beta12, 0.0, 1.0E3, log=T)
  beta11.priors <- dnormv(beta11, 0.0, 1.0E3, log=T)
  beta22.priors <- dnormv(beta22, 0.0, 1.0E3, log=T)
#
  LP <- LL + beta0.priors + #beta1.priors + 
        beta2.priors + beta12.priors + beta11.priors + beta22.priors
#
 
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rbinom(n=length(pi),size=n, prob=pi), parm=parm)
  return(Modelout)  
}
#
Model(parm, Data)
#
          r0 =c(0,0,1,3,0,0,4,0,0,2,0,0,0,0,1,0,0,1,0) 
          n0 =c(10,31,56,13,17,43,44,1,7,33,21,1,12,31,22,9,19,16,1)
          ST =c(1,1,1,1,1.7,1.7,1.7,1.7,2.2,2.2,2.2,2.2,2.8,2.8,2.8,4,4,4,4)
          HT =c(7,14,27,51,7,14,27,51,7,14,27,51,7,14,27,7,14,27,51)
# 
 parm <- c(beta0 =-3.9,# beta1 =-0.425,
           beta2 =0.11, beta12 =-0.09, beta11 =-0.45, beta22 =-0.005)
 Data <- list(N= length(r),mon.names=c("LogPosterior"), parm.names=names(parm),
              r0=r0, n0=n0, ST=ST, HT=HT)
# Run LaplacesDemon
out1 <- LaplaceApproximation(Model, Data=Data, parm=parm, Iterations=1000, Method="HAR")
#
parm <- as.data.frame(out1$Summary2)$Mean[1:5]
out2 <- LaplacesDemon.hpc(Model, Data=Data, Initial.Values=parm, Chains=3, CPUs=2,
                     Iterations=500000, Status=2000, Thinning=50)
#
# 
parm <- as.data.frame(out2[[3]]$Summary2)$Mean[1:5]
out <- LaplacesDemon(Model, Data=Data, Initial.Values=parm,
                         Iterations=100000, Status=2000, Thinning=50)
#
as.data.frame(out$Summary2)$Mean[1:5]
plot(out, BurnIn=50, Data, PDF=F)
######################################################################################
# read in joint posterior draws of the regression coefficients:
# beta0, beta1, beta2, beta12, beta11, beta22
# joint posterior draws in file "posteriorDraws.txt" from WinBUGS
# 100,000 draws for each parameter
#############################################
#out$Posterior1
data = out$Posterior1
#
beta0 = data[,1]
#beta1 = data[,2]
beta2 = data[,2]
beta12 = data[,3]
beta11 = data[,4]
beta22 = data[,5]
ndraws = length(beta0)
# separate function for calculating Bayesian chi-squared test statistics
ChiSquareTestStats <- function(N,K,x){
  Emk = N/K
  a = seq(0, 1, length.out=(K+1))
  mk = table(cut(x, breaks=a,include.lowest=T))
  return(sum((mk-Emk)^2/Emk))        }  # K = number of bins
# --------------------------------------------#
# p.value = p-value posterior distribution draws
K = floor(N^0.4)
testStat = rep(0,ndraws)
for(i in 1:ndraws){
  pi0 = exp(beta0[i]+beta1[i]*X1+beta2[i]*X2+beta12[i]*X1*X2+beta11[i]*X1^2+
           beta22[i]*X2^2)
  pi = pi0/(1+pi0)
  u = runif(N, pbinom(r-1, n, pi), pbinom(r, n, pi))
  testStat[i] = ChiSquareTestStats(N,K,u)
                   }
p.value = 1 - pchisq(testStat, (K-1))
median.p.value = median(p.value)
print(c("median posterior p-value ", median.p.value))
###############################################################################
# d = dimension of vector of regression coefficients
# BIC = Bayesian Information Criterion
d = 6
# compute regression coefficient posterior means
beta0.mean = mean(beta0)
beta1.mean = mean(beta1)
beta2.mean = mean(beta2)
beta12.mean = mean(beta12)
beta11.mean = mean(beta11)
beta22.mean = mean(beta22)
exp0 = exp(beta0.mean+beta1.mean*X1+beta2.mean*X2+beta12.mean*X1*X2+
      beta11.mean*X1^2+beta22.mean*X2^2)
pi.mean = exp0/(1+exp0)
BIC = -2*log(prod(dbinom(r,n,pi.mean)))+log(N)*d
print(c("BIC ", BIC))
###################################################################################
# R Code for Contour Plots
# separate function for finding matrix indices
matrix.indices <- function(a,b){
  #a is the number that the indices is wanted for
  #b is the matrix to extract the indices from
for(i in 1:nrow(b)){ #i is for rows
  for(j in 1:ncol(b)){ #j is for columns
    if(b[i,j] == a) return(c(i,j))
  }
}
}
  # ------------------------------------------------
  # cont.m = matrix of posterior median values
  # cont.q = matrix for 90th percentiles
  # x.m = matrix of indices for min in median
  # x.q = matrix of indices for min in 90th percentiles
  L1 = seq(1,4,.1)
  L1.c = L1 - mean(ST)
  L2 = seq(7,51,.5)
  L2.c = L2 - mean(HT)
  L1.n = length(L1)
  L2.n = length(L2)
  cont.m = matrix(0,nrow = L1.n, ncol = L2.n)
  cont.q = matrix(0,nrow = L1.n, ncol = L2.n)
  for(i in 1:L1.n){
    for(j in 1:L2.n){
      exp0 = exp(beta0+beta1*L1.c[i]+beta2*L2.c[j]+beta12*L1.c[i]*L2.c[j]+
            beta11*L1.c[i]^2+beta22*L2.c[j]^2)
      pi = exp0/(1+exp0)
      cont.m[i,j] = median(pi)
      cont.q[i,j] = quantile(pi, probs = .90)}}
  x.m = matrix.indices(min(cont.m), cont.m)
  x.q = matrix.indices(min(cont.q), cont.q)
  y = cbind(L1, 1:nrow(cont.m))
  z = cbind(L2, 1:ncol(cont.m))
  contour(L1,L2,cont.m,levels=seq(0,0.2,.005),main = "",xlab = "Soaking
Time",ylab="Heating Time")
  points(y[x.m[1],1],z[x.m[2],1],pch="*",cex=2)
  
  contour(L1,L2,cont.q,levels = seq(0,1,.005),main = "",xlab = "Soaking
Time",ylab = "Heating Time")
  points (y[x.q[1],1],z[x.q[2] ,1] ,pch = "*", cex = 2)
#  R Code for Prediction
  # predictive distribution at (STstar, HTstar) for predn test units
  # STstar = soaking time factor level
  # HTstar = heating time factor level
  #p redn = number of test units
  # predx = draws from predictive distribution of X for
  # X ~ Binomial(predn, pi),
  # where X is the number of failed units and
  # pi = posterior draws of failure probability at (STstar, HTstar)

  predn  = 100
  STstar = 1
  HTstar = 7
  X1star = STstar - mean(ST)
  X2star = HTstar - mean(HT)
  predn.vec = rep(predn, ndraws)
  exp0 = exp(beta0+beta1*X1star+beta2*X2star+beta12*X1star*X2star+beta11*X1star^2+
        beta22*X2star^2)
  pi = exp0 / (1+exp0)
  predx = rbinom(ndraws,predn.vec,pi)
  # table of predictive values and predictive probabilities
  print("predictive values and probabilities")
  table(predx)/length(predx)
  ##########################################################################
  
  
  
  













