###########################################################################
# LaplaceApproximation                                                    #
#                                                                         #
# The purpose of the LaplaceApproximation function is to maximize the     #
# logarithm of the unnormalized joint posterior distribution of a         #
# Bayesian model with a deterministic approximation that is a resilient   #
# backpropagation algorithm (among other selectable algorithms) to        #
# estimate the posterior modes and variances.                             #
###########################################################################

LaplaceApproximation <- function(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="LBFGS", Samples=1000, sir=TRUE,
     Stop.Tolerance=1.0E-5)
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     LA.call <- match.call()
     if(missing(Model)) stop("Model is a required argument.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("Data is a required argument.")
     if(missing(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to LaplaceApproximation().\n")
          parm <- rep(0, length(Data$parm.names))}
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     if({Interval <= 0} | {Interval > 1}) Interval <- 1.0E-6
     Iterations <- round(Iterations)
     if(Iterations < 10) {Iterations <- 10}
     else if(Iterations > 1000000) {Iterations <- 1000000}
     if({Method != "AGA"} & {Method != "HAR"} & {Method != "LBFGS"} &
          {Method != "NM"} & {Method != "Rprop"} & {Method != "SOMA"})
          stop("Method is unknown.")
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out,1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of apply functions\n")
          cat(     "were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are 'vectorized'.\n")
     }
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are 'vectorized'.\n")
     }
     ### Sample Size of Data
     if(!is.null(Data$n)) if(length(Data$n) == 1) N <- Data$n
     if(!is.null(Data$N)) if(length(Data$N) == 1) N <- Data$N
     if(!is.null(Data$y)) N <- nrow(matrix(Data$y))
     if(!is.null(Data$Y)) N <- nrow(matrix(Data$Y))
     if(!is.null(N)) cat("Sample Size: ", N, "\n")
     else stop("Sample size of Data not found in n, N, y, or Y.")
     ###########################  Preparation  ############################
     ### Attempt prior to starting...
     m.old <- Model(parm, Data)
     if(!is.list(m.old)) stop("Model must return a list.")
     if(length(m.old) != 5) stop("Model must return five components.")
     if(any(names(m.old) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.")
     if(length(m.old[["LP"]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(parm), length(m.old[["parm"]])))
          stop("The number of initial values and parameters differs.")
     if(!is.finite(m.old[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          if(!is.null(Data$PGF))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data)
          m.old <- Model(Initial.Values, Data)
          }
     if(!is.finite(m.old[["LP"]]))
          stop("The posterior is non-finite.")
     if(!is.finite(m.old[["Dev"]]))
          stop("The deviance is non-finite.")
     parm <- m.old[["parm"]]
     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation begins...\n")
     if(Method == "AGA") {
          LA <- AGA(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)}
     if(Method == "HAR") {
          LA <- HAR(Model, parm, Data, Iterations, Stop.Tolerance, m.old)}
     if(Method == "LBFGS") {
          LA <- LBFGS(Model, parm, Data, Iterations, Stop.Tolerance, m.old)}
     if(Method == "NM") {
          LA <- NM(Model, parm, Data, Iterations, Stop.Tolerance, m.old)}
     if(Method == "Rprop") {
          LA <- Rprop(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old)}
     if(Method == "SOMA") {
          LA <- SOMA(Model, parm, Data,
               options=list(maxiter=Iterations, tol=Stop.Tolerance))}
     Dev <- as.vector(LA$Dev)
     iter <- LA$iter
     parm.len <- LA$parm.len
     parm.new <- LA$parm.new
     parm.old <- LA$parm.old
     post <- LA$post
     Step.Size <- LA$Step.Size
     tol.new <- LA$tol.new
     if(iter == 1) stop("LaplaceApproximation stopped at iteration 1.")
     converged <- ifelse(tol.new <= Stop.Tolerance, TRUE, FALSE)
     ### Column names to samples
     if(ncol(post) == length(Data$parm.names))
          colnames(post) <- Data$parm.names
     rownames(post) <- 1:{nrow(post)}
     ########################  Covariance Matirx  #########################
     Approx.Hessian <- Hessian(Model, parm.new, Data)
     Inverse.test <- try(-as.inverse(Approx.Hessian), silent=TRUE)
     if(!inherits(Inverse.test, "try-error")) {
          VarCov <- Inverse.test
          diag(VarCov) <- ifelse(diag(VarCov) <= 0,
               .Machine$double.eps, diag(VarCov))}
     else {
          cat("\nWARNING: Failure to solve matrix inversion of ",
               "Approx. Hessian.\n", sep="")
          cat("NOTE: Identity matrix is supplied instead.\n")
               VarCov <- diag(parm.len)}
     #################  Sampling Importance Resampling  ##################
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Sampling from Posterior with Sampling Importance Resampling\n")
          posterior <- SIR(Model, Data, mu=parm.new, Sigma=VarCov,
               n=Samples)
          Mon <- matrix(0, nrow(posterior), length(Data$mon.names))
          dev <- rep(0, nrow(posterior))
          for (i in 1:nrow(posterior)) {
               mod <- Model(posterior[i,], Data)
               dev[i] <- mod[["Dev"]]
               Mon[i,] <- mod[["Monitor"]]
               }
          colnames(Mon) <- Data$mon.names}
     else {
          if({sir == TRUE} & {converged == FALSE})
               cat("Posterior samples are not drawn due to Converge=FALSE\n")
          posterior <- NA; Mon <- NA}
     #####################  Summary, Point-Estimate  ######################
     cat("Creating Summary from Point-Estimates\n")
     Summ1 <- matrix(NA, parm.len, 4, dimnames=list(Data$parm.names,
          c("Mode","SD","LB","UB")))
     Summ1[,1] <- parm.new
     Summ1[,2] <- sqrt(diag(VarCov))
     Summ1[,3] <- parm.new - 2*Summ1[,2]
     Summ1[,4] <- parm.new + 2*Summ1[,2]
     ###################  Summary, Posterior Samples  ####################
     Summ2 <- NA
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Creating Summary from Posterior Samples\n")
          Summ2 <- matrix(NA, ncol(posterior), 7,
               dimnames=list(Data$parm.names,
                    c("Mean","SD","MCSE","ESS","LB","Median","UB")))
          Summ2[,1] <- colMeans(posterior)
          Summ2[,2] <- apply(posterior, 2, sd)
          Summ2[,3] <- Summ2[,2] / sqrt(nrow(posterior))
          Summ2[,4] <- rep(nrow(posterior), ncol(posterior))
          Summ2[,5] <- apply(posterior, 2, quantile, c(0.025))
          Summ2[,6] <- apply(posterior, 2, quantile, c(0.500))
          Summ2[,7] <- apply(posterior, 2, quantile, c(0.975))
          Deviance <- rep(0, 7)
          Deviance[1] <- mean(dev)
          Deviance[2] <- sd(dev)
          Deviance[3] <- sd(dev) / sqrt(nrow(posterior))
          Deviance[4] <- nrow(posterior)
          Deviance[5] <- as.numeric(quantile(dev, probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(dev, probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(dev, probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:ncol(Mon)) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[,j])
               Monitor[2] <- sd(as.vector(Mon[,j]))
               Monitor[3] <- sd(as.vector(Mon[,j])) / sqrt(nrow(Mon))
               Monitor[4] <- nrow(Mon)
               Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
                    na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]
               }
          }
     ###############  Logarithm of the Marginal Likelihood  ###############
     LML <- LML(Model, Data, parm.new, method="LME2")
     LML <- list(LML=NA, VarCov=VarCov)
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          lml <- LML(theta=posterior, LL=(dev*(-1/2)), method="NSIS")
          LML[[1]] <- lml[[1]]}
     else if({sir == FALSE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(Model, Data, Modes, method="LME2")}
     colnames(VarCov) <- rownames(VarCov) <- Data$parm.names
     time2 <- proc.time()
     #############################  Output  ##############################
     LA <- list(Call=LA.call,
          Converged=converged,
          Covar=VarCov,
          Deviance=as.vector(Dev),
          History=post,
          Initial.Values=parm,
          Iterations=iter,
          LML=LML[[1]],
          LP.Final=as.vector(Model(parm.new, Data)[["LP"]]),
          LP.Initial=Model(parm, Data)[["LP"]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60, 2),
          Monitor=Mon,
          Posterior=posterior,
          Step.Size.Final=Step.Size,
          Step.Size.Initial=1,
          Summary1=Summ1,
          Summary2=Summ2,
          Tolerance.Final=tol.new,
          Tolerance.Stop=Stop.Tolerance)
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
     }
AGA <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old)
     {
     alpha.star <- 0.234
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     parm.new <- parm.old <- m.old[["parm"]]
     names(parm.new) <- names(parm.old) <- Data$parm.names
     tol.new <- 1
     Step.Size  <- 1 / parm.len
     post <- matrix(parm.new, 1, parm.len)
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Approximate Truncated Gradient
          approx.grad <- partial(Model, parm.old, Data, Interval)
          approx.grad <- interval(approx.grad, -1000, 1000, reflect=FALSE)
          ### Proposal
          parm.new <- parm.old + Step.Size * approx.grad
          if(any(!is.finite(parm.new))) parm.new <- parm.old
          m.new <- Model(parm.new, Data)
          tol.new <- sqrt(sum({parm.new - parm.old}^2))
          if(!is.finite(m.new[["LP"]])) {
               m.new <- m.old
               parm.new <- parm.old}
          ### Accept/Reject and Adapt
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               parm.new <- m.new[["parm"]]
               Step.Size <- Step.Size + (Step.Size / (alpha.star *
                    (1 - alpha.star))) * (1 - alpha.star) / iter
               }
          else {
               m.new <- m.old
               parm.new <- parm.old
               Step.Size <- abs(Step.Size - (Step.Size / (alpha.star *
                    (1 - alpha.star))) * alpha.star / iter)
               }
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
HAR <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old)
     {
     alpha.star <- 0.234
     tau <- 1
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     parm.new <- parm.old <- parm
     names(parm.new) <- names(parm.old) <- Data$parm.names
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Propose new parameter values
          theta <- rnorm(parm.len)
          d <- theta / sqrt(sum(theta*theta))
          Step.Size <- runif(1,0,tau)
          parm.new <- parm.old + Step.Size * d
          m.new <- Model(parm.new, Data)
          tol.new <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
          if(!is.finite(m.new[["LP"]])) m.new <- m.old
          if(!is.finite(m.new[["Dev"]])) m.new <- m.old
          if(any(!is.finite(m.new[["Monitor"]]))) m.new <- m.old
          ### Accept/Reject and Adapt
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               parm.new <- m.new[["parm"]]
               tau <- tau + (tau / (alpha.star *
                    (1 - alpha.star))) * (1 - alpha.star) / iter
               }
          else {
               m.new <- m.old
               parm.new <- parm.old
               tau <- abs(tau - (tau / (alpha.star *
                    (1 - alpha.star))) * alpha.star / iter)
               }
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
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
          return(out)
          }
     for (iter in 1:Iterations) {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### LBFGS
          Fit <- optim(par=parm.new, fn=ModelWrapper,
               method="L-BFGS-B", control=list(fnscale=-1, maxit=1))
          m.new <- Model(Fit$par, Data)
          tol <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
          if(!is.finite(m.new[["LP"]]) | {m.new[["LP"]] < m.old[["LP"]]})
               m.new <- m.old
          parm.new <- m.new[["parm"]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=0,
          tol.new=tol)
     return(LA)
     }
NM <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
     ...)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     n <- length(as.vector(parm))
     parm.old <- m.old[["parm"]]
     Step.Size <- 1
     if(!is.numeric(parm) || length(parm) < 2)
          stop("Parameters must be a numeric vector of length > 1.")
     # simplex vertices around parm
     V <- t(1/(2*n) * cbind(diag(n), rep(-1, n)) + parm)
     P <- Q <- c()
     # Function values at vertices
     Y <- numeric(n+1)
     for (j in 1:(n+1)) {
          Mo <- Model(V[j, ], Data)
          Y[j] <- Mo[["LP"]] * -1
          V[j,] <- Mo[["parm"]]
          }
     ho <- lo <- which.min(Y)
     li <- hi <- which.max(Y)
     for (j in 1:(n+1)) {
          if(j != lo && j != hi && Y[j] <= Y[li]) li <- j
          if(j != hi && j != lo && Y[j] >= Y[ho]) ho <- j}
     iter <- 0
     while(Y[hi] > Y[lo] + Stop.Tolerance && iter < Iterations) {
          S <- numeric(n)
          for (j in 1:(n+1)) S <- S + V[j,1:n]
          M <- (S - V[hi,1:n])/n
          R <- 2*M - V[hi,1:n]
          Mo <- Model(R, Data)
          yR <- Mo[["LP"]] * -1
          R <- Mo[["parm"]]
          if(yR < Y[ho]) {
               if(Y[li] < yR) {
                    V[hi,1:n] <- R
                    Y[hi] <- yR
                    }
               else {
                    E <- 2*R - M
                    Mo <- Model(E, Data)
                    yE <- Mo[["LP"]] * -1
                    E <- Mo[["parm"]]
                    if(yE < Y[li]) {
                         V[hi,1:n] <- E
                         Y[hi] <- yE
                         }
                    else {
                         V[hi,1:n] <- R
                         Y[hi] <- yR
                         }
                    }
               }
          else {
               if(yR < Y[hi]) {
                    V[hi,1:n] <- R
                    Y[hi] <- yR}
               C <- (V[hi,1:n] + M) / 2
               Mo <- Model(C, Data)
               yC <- Mo[["LP"]] * -1
               C <- Mo[["parm"]]
               C2 <- (M + R) / 2
               Mo <- Model(C2, Data)
               yC2 <- Mo[["LP"]] * -1
               C2 <- Mo[["parm"]]
               if(yC2 < yC) {
                    C <- C2
                    yC <- yC2}
               if(yC < Y[hi]) {
                    V[hi,1:n] <- C
                    Y[hi] <- yC
                    }
               else {
                    for (j in 1:(n+1)) {
                         if(j != lo) {
                              V[j,1:n] <- (V[j,1:n] + V[lo,1:n]) / 2
                              Z <- V[j,1:n]
                              Mo <- Model(Z, Data)
                              Y[j] <- Mo[["LP"]] * -1
                              Z <- Mo[["parm"]]
                              V[j,1:n] <- Z}
                         }
                    }
               }
          ho <- lo <- which.min(Y)
          li <- hi <- which.max(Y)
          for (j in 1:(n+1)) {
               if(j != lo && j != hi && Y[j] <= Y[li]) li <- j
               if(j != hi && j != lo && Y[j] >= Y[ho]) ho <- j
               }
          iter <- iter + 1
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          P <- rbind(P, V[lo, ])
          Q <- c(Q, Y[lo])
          Dev <- rbind(Dev, Model(V[lo,], Data)[["Dev"]])
          }
     snorm <- 0
     for (j in 1:(n+1)) {
          s <- abs(V[j] - V[lo])
          if(s >= snorm) snorm <- s}
     V0 <- V[lo, 1:n]
     y0 <- Y[lo]
     dV <- snorm
     dy <- abs(Y[hi] - Y[lo])
     Dev <- Dev[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=V0,
          parm.old=parm.old, post=P, Step.Size=Step.Size,
          tol.new=Y[hi] - Y[lo])
     return(LA)
     }
Rprop <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old)
     {
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.len <- length(as.vector(parm))
     approx.grad.old <- approx.grad.new <- rep(0, parm.len)
     parm.old <- m.old[["parm"]]
     parm.new <- parm.old - 0.1 #First step
     names(parm.new) <- names(parm.old) <- Data$parm.names
     tol.new <- 1
     post <- matrix(parm.new, 1, parm.len)
     Step.Size <- rep(0.0125, parm.len)
     for (iter in 1:Iterations) {
          approx.grad.old <- approx.grad.new
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Approximate Truncated Gradient
          approx.grad.new <- partial(Model, parm.old, Data, Interval)
          approx.grad.new <- interval(approx.grad.new, -1000, 1000,
               reflect=FALSE)
          ### Determine if Gradients Changed Sign
          change <- (approx.grad.old >= 0) == (approx.grad.new >= 0)
          ### Adjust weight (step size) based on sign change
          Step.Size[change == TRUE] <- Step.Size[change == TRUE] * 0.5
          Step.Size[change == FALSE] <- Step.Size[change == FALSE] * 1.2
          Step.Size <- interval(Step.Size, 0.0001, 50, reflect=FALSE)
          ### Propose new state based on weighted approximate gradient
          parm.new <- parm.old + Step.Size * approx.grad.new
          m.new <- Model(parm.new, Data)
          tol.new <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
          if(!is.finite(m.new[["LP"]]) | {m.new[["LP"]] < m.old[["LP"]]}) {
               p.order <- sample(1:length(parm.new))
               parm.temp <- parm.old
               for (i in 1:length(p.order)) {
                    parm.temp[p.order[i]] <- parm.new[p.order[i]]
                    m.new <- Model(parm.temp, Data)
                    if(m.new[["LP"]] < m.old[["LP"]])
                         parm.temp[p.order[i]] <- parm.old[p.order[i]]}
               m.new <- Model(parm.temp, Data)}
          parm.new <- m.new[["parm"]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[["Dev"]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
SOMA <- function (Model, parm, Data, bounds, options=list())
     {
     ### Initial Checks
     if(missing(bounds)) {
          bounds <- list(min=rep(-.Machine$double.xmax,
               length(parm)),
               max=rep(.Machine$double.xmax, length(parm)))
          }
     if(!(all(c("min","max") %in% names(bounds))))
          stop("Bounds list must contain \"min\" and \"max\" vector elements")
     if(length(bounds$min) != length(bounds$max))
          stop("Bounds are not of equal length.")
    ### Option Defaults
    defaultOptions <- list(pathLength=3, stepLength=0.11,
         perturbationChance=0.1, tol=1e-3, maxiter=1000,
         populationSize=10)
    defaultsNeeded <- setdiff(names(defaultOptions), names(options))
    spuriousOptions <- setdiff(names(options), names(defaultOptions))
    options[defaultsNeeded] <- defaultOptions[defaultsNeeded]
    if(length(spuriousOptions) > 0)
          warning("The following specified options are unused: ",
               paste(spuriousOptions,collapse=", "))
    ### Setup
    nParams <- length(bounds$min)
    nParamsTotal <- nParams * options$populationSize
    steps <- seq(0, options$pathLength, options$stepLength)
    nSteps <- length(steps)
    steps <- rep(steps, each=nParamsTotal)
    post <- matrix(parm, 1, nParams)
    ### Create Population
    population <- matrix(parm, nrow=nParams, ncol=options$populationSize)
    for (i in 2:ncol(population)) {
         if(!is.null(Data$PGF))
              population[,i] <- GIV(Model, Data, PGF=TRUE)
         else population[,i] <- GIV(Model, Data)}
    ### Calculate initial LP and Dev per individual in population
    temp <- apply(population, 2, function(x) Model(x, Data))
    population <- sapply(temp, function(l) l$parm)
    LP <- sapply(temp, function(l) l$LP)
    Dev <- sapply(temp, function(l) l$Dev)
    iteration <- 0
    DevHistory <- LPHistory <- numeric(0)
    if((max(LP) - min(LP)) < options$tol) {
         population <- matrix(runif(nParamsTotal,-10,10), nrow=nParams,
              ncol=options$populationSize)}
    if(all(!is.finite(LP))) stop("All individuals have non-finite LPs.")
    ### Evolution Begins
    repeat {
         ### Find the current leader
         leaderIndex <- which.max(LP)
         leaderDev <- Dev[leaderIndex]
         leaderLP <- LP[leaderIndex]
         ### Check termination criteria
         if(iteration == options$maxiter) break
         LPdiff <- max(LP) - min(LP)
         if(!is.finite(LPdiff)) LPdiff <- 1
         if(LPdiff < options$tol) break
         DevHistory <- c(DevHistory, leaderDev)
         LPHistory <- c(LPHistory, leaderLP)
         ### Find the migration direction for each individual
         directionsFromLeader <- apply(population, 2, "-",
              population[,leaderIndex])
         ### Establish which parameters will be changed
         toPerturb <- runif(nParamsTotal) < options$perturbationChance
         # Second line here has a minus, so directions are away from leader
         populationSteps <- array(rep(population, nSteps),
              dim=c(nParams, options$populationSize, nSteps))
         populationSteps <- populationSteps -
              steps * rep(directionsFromLeader * toPerturb, nSteps)
         ### Replace out-of-bounds parameters with random valid values
         outOfBounds <- which(populationSteps < bounds$min |
              populationSteps > bounds$max)
         randomSteps <- array(runif(nParamsTotal*nSteps),
              dim=c(nParams, options$populationSize, nSteps))
         randomSteps <- randomSteps * (bounds$max-bounds$min) + bounds$min
         populationSteps[outOfBounds] <- randomSteps[outOfBounds]
         ### Values over potential locations
         temp <- apply(populationSteps, 2:3, function(x) Model(x, Data))
         population <- sapply(temp, function(l) l$parm)
         LP <- sapply(temp, function(l) l$LP)
         LP <- matrix(LP, length(LP) / nSteps, nSteps)
         Dev <- sapply(temp, function(l) l$Dev)
         Dev <- matrix(Dev, length(Dev) / nSteps, nSteps)
         individualBestLocs <- apply(LP, 1, which.max)
         ### Migrate each individual to its best new location, and update LP
         indexingMatrix <- cbind(seq_len(options$populationSize),
              individualBestLocs)
         population <- t(apply(populationSteps, 1, "[", indexingMatrix))
         LP <- LP[indexingMatrix]
         Dev <- Dev[indexingMatrix]
         post <- rbind(post, population[,leaderIndex])
         iteration <- iteration + 1
         if(iteration %% round(options$maxiter / 10) == 0)
              cat("Iteration: ", iteration, " of ", options$maxiter, "\n")
         }
    ### Output
    LA <- list(Dev=DevHistory, 
         iter=iteration,
         parm.len=nParams,
         parm.new=population[,leaderIndex],
         parm.old=parm,
         post=post[-1,],
         Step.Size=options$stepLength,
         tol.new=signif(LPdiff,3))
    return(LA)
    }
#End
