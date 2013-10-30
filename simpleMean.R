# Simple normal mean model in LaplacesDemon
# Generate two samples of body mass measurements of male peregrines
y1000 <- rnorm(n = 1000, mean = 100, sd = 10)  # Sample of 1000 birds
###==========================================================
mean(y1000)
###==========================================================
lm0 <- lm(y1000~1)
sd(y1000)
summary(lm0)
###=========================================================

population.sd <- 1 

for(i in 1:10){
mu <- 1:3000
la00 <- sapply(mu,function(xx)sum(dnorm(y1000, xx, population.sd, log=TRUE)))
mu <- mu[which.max(la00)]
population.sd <- 1:100
d01 <- sapply(population.sd,function(xx)sum(dnorm(y1000, mu, xx, log=TRUE)))
population.sd <- population.sd[which.max(d01)]
}
c(mean=mu,sd=population.sd)

plot(la00)
plot(d01)
 
###===============================================================
### Random walk MCMC for binomial proportion
############################################
# Parameters
parm <- c(1,10)
population.mean <- parm[1]
population.sd <- parm[2]
 
# Prior density
population.mean.prior <- dunif(population.mean, 0, 5000)
population.sd.prior <- dunif(population.sd, 0, 100)

# Log-Likelihood
LL <- sum(dnorm(y1000, population.mean, population.sd, log=TRUE))

# Log-Posterior
LP <- LL + population.mean.prior + population.sd.prior
Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP), 
                 yhat=rnorm(length(y1000), population.mean, population.sd),
                 parm=c(rnorm(1,600), runif(1, 1, 30))
                 )

###===============================================================

# Load library 
library(LaplacesDemon)
# Model specification
Model <- function(parm, Data)
{
  # Parameters
  population.mean <- parm[1]
  population.sd <- parm[2]
  # Prior density
  population.mean.prior <- dunif(population.mean, 0, 5000)
  population.sd.prior <- dunif(population.sd, 0, 100)
  # Log-Likelihood
  mu <- population.mean
  LL <- sum(dnorm(Data$mass, mu, population.sd, log=TRUE))
  # Log-Posterior
  LP <- LL + population.mean.prior + population.sd.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP), yhat=rnorm(Data$N, mu, population.sd), parm=parm)
  return(Modelout)
}
# Prepare the data 
parm.names <- c("population.mean", "population.sd")
Data <- list(mass=y1000, N=length(y1000), mon.names=c("LP"), parm.names=parm.names)
# Initial values
Initial.Values <- c(
  rnorm(1,600),   # population.mean
  runif(1, 1, 30) # population.sd
)
# MCMC settings
ni <- 50000  	# Number of draws from posterior (for each chain)
st <- 1000 		# Steps when status message should be given
nt <- 50 	# Thinning rate
# Run LaplacesDemon
out <- LaplacesDemon(Model, Data=Data, Initial.Values, Iterations=ni, Status=st, Thinning=nt)
# Have a look at some summary statistics
out
# Plotting output
plot(out, BurnIn=100, Data, PDF=T, Parms=c("population.mean", "population.sd"))