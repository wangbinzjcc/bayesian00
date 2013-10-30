###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################

LaplacesDemon <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=100000, Status=1000, Thinning=100, Algorithm="RWM",
     Specs=NULL, ...)
     
   ############################  Begin MCMC  ############################
    if(Algorithm == "Adaptive Hamiltonian Monte Carlo") {
          mcmc.out <- AHMC( ... )}
     else stop("The algorithm is unrecognized.")
     #########################  MCMC is Finished  #########################
     Acceptance <- mcmc.out$Acceptance
     Dev <- mcmc.out$Dev
     DiagCovar <- mcmc.out$DiagCovar
     Mon <- mcmc.out$Mon
     thinned <- mcmc.out$thinned
     VarCov <- mcmc.out$VarCov
     remove(mcmc.out)
     colnames(DiagCovar) <- Data$parm.names
     thinned <- matrix(thinned[-1,], nrow(thinned)-1, ncol(thinned))
     Dev <- matrix(Dev[-1,], nrow(Dev)-1, 1)
     Mon <- matrix(Mon[-1,], nrow(Mon)-1, ncol(Mon))
     if(is.matrix(VarCov) & !is.list(VarCov)) {
          colnames(VarCov) <- rownames(VarCov) <- Data$parm.names}
     else if(is.vector(VarCov) & !is.list(VarCov)) {
          names(VarCov) <- Data$parm.names}
     thinned.rows <- nrow(thinned)
     ### Warnings (After Updating)
     if(any(Acceptance == 0))
          cat("\nWARNING: All proposals were rejected.\n")
     ### Real Values
     thinned <- ifelse(!is.finite(thinned), 0, thinned)
     Dev <- ifelse(!is.finite(Dev), 0, Dev)
     Mon <- ifelse(!is.finite(Mon), 0, Mon)
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n")
     if(thinned.rows %% 10 == 0) thinned2 <- thinned
     if(thinned.rows %% 10 != 0) thinned2 <- thinned[1:(10*trunc(thinned.rows/10)),]
     HD <- BMK.Diagnostic(thinned2, batches=10)
     Ind <- 1 * (HD > 0.5)
     BurnIn <- thinned.rows
     batch.list <- seq(from=1, to=nrow(thinned2), by=floor(nrow(thinned2)/10))
     for (i in 1:9) {
          if(sum(Ind[,i:9]) == 0) {
               BurnIn <- batch.list[i] - 1
               break
               }
          }
     Stat.at <- BurnIn + 1
     rm(batch.list, HD, Ind, thinned2)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n")
     acf.temp <- matrix(1, trunc(10*log10(thinned.rows)), LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=nrow(acf.temp), plot=FALSE)
          acf.temp[,j] <- abs(temp0$acf[2:{nrow(acf.temp)+1},,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning}
     Rec.Thin <- ifelse(is.na(Rec.Thin), nrow(acf.temp), Rec.Thin)
     ### Assess ESS for all deviance and monitor samples
     ESS2 <- ESS(Dev)
     ESS3 <- ESS(Mon)
     ### Assess ESS for stationary samples
     if(Stat.at < thinned.rows) {
          ESS4 <- ESS(thinned[Stat.at:thinned.rows,])
          ESS5 <- ESS(Dev[Stat.at:thinned.rows,])
          ESS6 <- ESS(Mon[Stat.at:thinned.rows,])}
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- apply(thinned, 2, sd)
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)), silent=TRUE)
     if(inherits(temp, "try-error"))
          temp <- MCSE(as.vector(Dev), method="sample.variance")
     Deviance[3] <- temp
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
          if(inherits(temp, "try-error")) 
               temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.5,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data$mon.names[j]}
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(Stat.at < thinned.rows) {
          thinned2 <- matrix(thinned[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(thinned))
          Dev2 <- matrix(Dev[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Dev))
          Mon2 <- matrix(Mon[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- apply(thinned2, 2, sd)
          Summ2[,3] <- 0
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=TRUE)
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else Summ2[i,3] <- MCSE(thinned2[,i],
                    method="sample.variance")}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(as.vector(Dev2), method="sample.variance")
          Deviance[3] <- temp
          Deviance[4] <- ESS5
          Deviance[5] <- as.numeric(quantile(Dev2, probs=0.025,
               na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev2, probs=0.500,
               na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev2, probs=0.975,
               na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon2[,j])
               Monitor[2] <- sd(as.vector(Mon2[,j]))
               temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
               if(inherits(temp, "try-error"))
                    temp <- MCSE(as.vector(Mon[,j]),
                    method="sample.variance")
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]}
          }
     ### Column names to samples
     if(identical(ncol(Mon), length(Data$mon.names)))
          colnames(Mon) <- Data$mon.names
     if(identical(ncol(thinned), length(Data$parm.names))) {
          colnames(thinned) <- Data$parm.names}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(({Algorithm == "Affine-Invariant Ensemble Sampler"} |
          {Algorithm == "Componentwise Hit-And-Run Metropolis"} |
          {Algorithm == "Componentwise Slice"} |
          {Algorithm == "Delayed Rejection Metropolis"} |
          {Algorithm == "Hit-And-Run Metropolis"} | 
          {Algorithm == "Hamiltonian Monte Carlo"} |
          {Algorithm == "Independence Metropolis"} |
          {Algorithm == "Metropolis-within-Gibbs"} |
          {Algorithm == "No-U-Turn Sampler"} | 
          {Algorithm == "Random-Walk Metropolis"} |
          {Algorithm == "Reversible-Jump"} |
          {Algorithm == "Sequential Metropolis-within-Gibbs"} |
          {Algorithm == "Slice Sampler"} | 
          {Algorithm == "Tempered Hamiltonian Monte Carlo"} | 
          {Algorithm == "t-walk"}) & 
          {Stat.at < thinned.rows}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(theta=thinned2,
               LL=as.vector(Dev2)*(-1/2), method="NSIS")}
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance/Iterations,7),
          Adaptive=Adaptive,
          Algorithm=Algorithm,
          Call=LDcall,
          Covar=VarCov,
          CovarDHis=DiagCovar,
          Deviance=as.vector(Dev),
          DIC1=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          DIC2=if(Stat.at < thinned.rows) {
               c(mean(as.vector(Dev2)),
               var(as.vector(Dev2))/2,
               mean(as.vector(Dev2)) +
               var(as.vector(Dev2))/2)}
               else rep(NA,3),
          DR=DR,
          Initial.Values=Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60,2),
          Model=Model,
          Monitor=Mon,
          Parameters=LIV,
          Periodicity=Periodicity,
          Posterior1=thinned,
          Posterior2=if(Stat.at < thinned.rows) {
               thinned[Stat.at:thinned.rows,]}
               else thinned[thinned.rows,],
          Rec.BurnIn.Thinned=BurnIn,
          Rec.BurnIn.UnThinned=BurnIn*Thinning,
          Rec.Thinning=min(1000, max(Rec.Thin)),
          Status=Status,
          Summary1=Summ1,
          Summary2=Summ2,
          Thinned.Samples=thinned.rows,
          Thinning=Thinning)
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n")
     return(LaplacesDemon.out)
     } 

#########################################################

AM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, VarCov)
     {
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
          if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(post[1:iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }  
   
HMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, epsilon, L)
     {
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
  
 
RWM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, VarCov)
     {
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- Mo0[["parm"]]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
          if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
 
 
 
#End
