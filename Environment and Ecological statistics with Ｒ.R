########################################################
#   Environment and Ecological statistics with R
#                                 wangbinzjcc 2013-12-13
#########################################################

#########
#  95%  置信区间 confidence interval
###
n.sims <- 1000
n.size <- 30
inside <- 0
for(i in 1:n.sims){
  y <- rnorm(n.size, mean=2.05, sd=0.34)
  se <- sd(y)/sqrt(n.size)
  int.95 <- mean(y) + qt(c(.025, .975), n.size-1)*se
  inside<- inside + sum(int.95[1]<2.05 & int.95[2]>2.05)
  }

inside/n.sims

######
# bootstrapping  自举法
#
x <- c(94, 38, 23, 19, 99, 16,10,40,30,60, 141)
boot.sample <- sample(x, size=length(x), replace=T)
boot.mean <- mean(boot.sample)
boot.mean <- c()
B <- 2000
for(i in 1:B){
  boot.sample <- sample(x,size=length(x),T)
  boot.mean[i] <- mean(boot.sample)
  }
boot.se <- sd(boot.mean)
boot.se
se <- sd(x)/sqrt(length(x)); se

#######

set.seed(1)
t.sample <- rt(10000, df=23)
p.value <- mean(t.sample>2.34)
#
y <- log(rlnorm(25, 1.9,1))
n.sim <- 5000
n <- length(y)
y.bar <- mean(y)
s.hat <- sd(y)
theta.i <- s.hat * sqrt((n-1)/rchisq(n.sim,n-1))
mu.i <- rnorm(n.sim, y.bar, theta.i/sqrt(n))
y.tilde <- rnorm(n.sim, mu.i, theta.i)
Pr <- mean(y.tilde>log(10))
x.tilde <- exp(y.tilde)
mu.x <- mean(x.tilde)
sigma.x <- sd(x.tilde)
CI.x <- mu.x + qt(c(0.025, 0.975), df=25-1)* sigma.x/sqrt(25-1)
#####################
#  randomForest
###################




