
# Site occupancy model in LaplacesDemon
# Select sample sizes (spatial and temporal replication)
R <- 200 
T <- 3
# Determine process parameters
psi <- 0.8    # Occupancy probability
p <- 0.5      # Detection probability

# Create structure to contain counts
y <- matrix(NA, nrow = R, ncol = T)

# Ecological process: Sample true occurrence (z, yes/no) from a Bernoulli (occurrence probability = psi)
z <- rbinom(n = R, size = 1, prob = psi) # Latent occurrence state

# Observation process: Sample detection/nondetection observations from a Bernoulli(with p) if z=1
for (j in 1:T){
  y[,j] <- rbinom(n = R, size = 1, prob = z * p)
}

# Look at truth and at our imperfect observations
sum(z)                 # Realized occupancy among 200 surveyed sites
sum(apply(y, 1, max))  # Observed occupancy

# Load library 
library(LaplacesDemon)

# Model specification
Model <- function(parm, Data)
{
  # Parameters
  logit.psi <- parm[1]
  psi <- invlogit(logit.psi)
  logit.p <- parm[2]
  p <- invlogit(logit.p)
  
  # Prior density
  psi.prior <- dunif(psi, 0, 1, log=TRUE)
  p.prior <- dunif(p, 0, 1, log=TRUE)
  
  # Log-Likelihood
  LL.Z <- sum(dbern(z, psi, log=TRUE))
  LL1 <- sum(dbern(y[,1], z*p, log=TRUE))
  LL2 <- sum(dbern(y[,2], z*p, log=TRUE))
  LL3 <- sum(dbern(y[,3], z*p, log=TRUE))
  LL <- LL1 + LL2 + LL3 + LL.Z
  
  # Log-Posterior
  LP <- LL + psi.prior + p.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, psi, p), yhat=1, parm=parm)
  return(Modelout)
}

# Prepare the data 
parm.names <- c("logitpsi", "logitp")
Data <- list(y=y, mon.names=c("LP", "psi", "p"), parm.names=parm.names)

# Initial values
Initial.Values <- c(
  logit(runif(1,0,1)),   # psi
  logit(runif(1,0,1))    # p
  )

# Run LaplacesDemon
ni <- 10000  	# Number of draws from posterior (for each chain)
st <- 1000		# Steps when status message should be given
nt <- 1		# Thinning rate

out <- LaplacesDemon(Model, Data=Data, Initial.Values,
                     Iterations=ni, Status=st, Thinning=nt)

# Have a look at some summary statistics
out
Consort(out)

out <- LaplacesDemon(Model, Data=Data, Initial.Values,
                     Iterations=ni, Status=st, Thinning=nt,
                     Algorithm="AMM", Specs=list(Adaptive=11, Periodicity=220, w=0.05))

out <- LaplacesDemon(Model, Data=Data, Initial.Values,
                     Iterations=ni, Status=st, Thinning=nt,
                     Algorithm="NUTS", Specs=list(A=50, delta=0.6, epsilon=NULL))

# Have a look at some summary statistics
out
print(out$Summary1[c("psi", "p", "Deviance"),], 3)

# Plotting output
plot(out, BurnIn=8000, Data, PDF=TRUE, Parms=c("psi", "p"))

