#########################################################
#  Bayes LaplaceDemon  Mean    wangbinzjcc   2014-1-7
#########################################################
#getwd()
setwd('F://GitHub//bayesian00')
n <- 8000     # Number of females
mu1  <- 105   # Population mean for females
sigma1 <- 8.5 # Population SD for females
y <- rnorm(n, mean= mu1, sd=sigma1) # Data for females
hist(y)
plot(y)

#################################################################

mean(y)
sd(y)

boxplot(y)
#
tiff('histMeanBayes.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par(mar=c(4,4,1,1))
hist(y)
dev.off()
#
tiff('PlotMeanBayes.tiff',
     width = 3000, height = 3000,res=600,compression = "lzw")
par(mar=c(4,4,1,1))
plot(y,xlab='模拟次数',cex=0.6)
abline(h=105, col=2, lwd=3)
dev.off()
##################################################################

require(LaplacesDemon)          #  R package
parm <-  c(mu1=105, log.sigm=log(8.5))    #  Initial.Values  
Data <- list(mon.names=c("sigmma"), parm.names=c("mean", "log.sigm"), y=y)

Model <- function(parm, Data){
  mu1 <- parm[1]; sigma1 <- exp(parm[2])
  LL <- sum(dnorm(x=Data$y, mean=mu1, sd=sigma1, log=T))
  mu1.prior <- dnormv(x=mu1, mean=0, var=1000, log=T)
  sigma.prior <- dgamma(sigma1, 25, log=T)
  LP <- LL + 2 * mu1.prior  + 2 * sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(sigma1),  
                   yhat= rnorm(length(Data$y), mu1, sigma1), parm=parm)
  return(Modelout)                   }

#out0 <- LaplacesDemon(Model, Data=Data, Initial.Values=parm, Iterations=10000)

ni <- 80000   # Number of draws from posterior (for each chain)
st <- 2000      # Steps when status message should be given
nt <- 50        # Thinning rate #  Abate autocorrelation


# Run LaplacesDemon
out <- LaplacesDemon(Model, Data=Data, Initial.Values=parm, 
                     Iterations=ni, Status=st, Thinning=nt)
# Have a look at some summary statistics
#
#
#
plot(out, BurnIn=10, Data, PDF=T)
#
#
mean(y)
sd(y)
##########################################################################
##########################################################################


setwd('F:\\GitHub\\DeadStandingTrees')
 
Xmat0 <- read.csv("Xmat.csv")

yDe <- read.csv("yDead.csv")
head(Xmat0)
head(yDe)
####################
require(LaplacesDemon)          #  R package
#

parm=c(alpha=1:3,beta=4:9, miu.a=10,miu.b=11,sigm.a_log=log(12),sigm.b_log=log(13))


parm=c(alpha=c(-0.05,-0.07,-0.1),beta=c(-0.06,-0.1,-0.06,-0.09,-0.1,-0.08), miu.a=0,miu.b=3.00,sigm.a_log=3,sigm.b_log=8)
parm.names <- as.parm.names(parm)
Data=list(x=Xmat0[,-1], parm.names=parm.names, mon.names=c("sigma.a","sigma.b"),
          y=yDe$x)

Model=function(parm, Data){
  # Priors
   alpha = parm[1:3]    # 3
   beta  = parm[4:9]    # 6
  #
   miu.a = parm[10]
   miu.b = parm[11]
    sigm.a = exp(parm[12])    # 1
     sigm.b = exp(parm[13])   # 1
  
  #
   alpha.prior <- sum(dnormv(alpha, miu.a, sigm.a, log=T))
   beta.prior <- sum(dnormv(beta, miu.b, sigm.b, log=T))
   miu.a.prior <- dnormv(miu.a, 0, 1000, log=T) 
   miu.b.prior <- dnormv(miu.b, 0, 1000, log=T) 
   sigm.a.prior <- dgamma(sigm.a, 25, log=T)
   sigm.b.prior <- dgamma(sigm.b, 25, log=T)  
  #
   mu= exp(rowSums(Data$x * c(alpha,beta)))
   LL <- sum(dpois(Data$y,mu,log=T))
  #
   LP <- LL + alpha.prior + beta.prior + miu.a.prior+ miu.b.prior + sigm.a.prior + sigm.b.prior
  
  # Model out
  yhat <- rpois(length(mu) ,mu)
  Modelout <- list(LP=LP, Dev= -2*LL, Monitor=c(sigm.a, sigm.b),
                 yhat=yhat, parm=parm)
  return(Modelout)              
}
###
out <- LaplacesDemon(Model, Initial.Values=parm, Data=Data,
                     Iterations=1000000, Status=200, Thinning=30)


#
# save.image("F:/GitHub/DeadStandingTrees/1000000timesBayes.RData")
# load("F:/GitHub/DeadStandingTrees/1000000timesBayes.RData")
out
# dput(out,"10000000timesLaplacesDenon")
# out0 <- dget("10000000timesLaplacesDenon")
summary(out)
rownames(summary(out))
#
o1 <- out[['Posterior1']]
o1 <- as.data.frame(o1)
head(o1)
aaa <- o1$alpha1
plot(aaa, type='o')
hist(aaa)
summary(aaa)
#
apply(o1, 2, quantile, probs=c(0.5, 0.05, 0.95))
out[['Summary2']]
###
quantile(aaa,probs=c(0.05,0.95))

out[[i]];i=i+1
summary(out[[16]] )
#
#
for(i in rownames(summary(out))){
  aa <- out[[i]]
try(write.csv(aa, paste(i, '.csv'))) 
}
#

out[['Summary2']]


####################################################################











