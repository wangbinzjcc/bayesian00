###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################
LaplacesDemon_wb <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=200 , Status=20 , Thinning=10, Algorithm="RWM",
     Specs=NULL, ...)     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     LDcall <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n")
     if(missing(Model)) stop("A function must be entered for Model.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data))
          stop("A list containing data must be entered for Data.")
     if(is.null(Data$mon.names)) stop("In Data, mon.names is NULL.")
     if(is.null(Data$parm.names)) stop("In Data, parm.names is NULL.")
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     if(missing(Initial.Values)) {
          cat("WARNING: Initial Values were not supplied.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(!identical(length(Initial.Values), length(Data$parm.names))) {
          cat("WARNING: The length of Initial Values differed from",
               "Data$parm.names.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     Iterations <- round(abs(Iterations))
     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="")}
     Status <- round(abs(Status))
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="")}
     Thinning <- round(abs(Thinning))
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="")}
     if(Algorithm %in% c("AHMC","AIES","AM","AMM","AMWG","CHARM","DEMC",
          "DRAM","DRM","Experimental","HARM","HMC","HMCDA","IM","INCA",
          "MWG","NUTS","RAM","RJ","RWM","SAMWG","Slice","SMWG","THMC",
          "twalk","USAMWG","USMWG")) {
          if(Algorithm == "AHMC") {
               Algorithm <- "Adaptive Hamiltonian Monte Carlo"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(all(Specs[["epsilon"]] == Specs[[1]])) {
                    epsilon <- as.vector(abs(Specs[[1]]))
                    if(length(epsilon) != length(Initial.Values))
                         cat("\nLength of epsilon is incorrect.\n")
                         epsilon <- rep(epsilon[1], length(Initial.Values))}
               if(Specs[["L"]] == Specs[[2]]) L <- abs(round(Specs[[2]]))
               if(L < 1) {
                    cat("\nL has been increased to its minimum: 1.\n")
                    L <- 1}
               Adaptive <- Iterations + 1
               DR <- 1
               if(Specs[["Periodicity"]] == Specs[[3]]) Periodicity <- Specs[[3]]}
          if(Algorithm == "AIES") {
               Algorithm <- "Affine-Invariant Ensemble Sampler"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 6) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1
               if(Specs[["Nc"]] == Specs[[1]]) Nc <- abs(round(Specs[[1]]))
               if(Nc < 3) Nc <- 3
               if(is.null(Specs[[2]]) | all(Specs[["Z"]] == Specs[[2]])) {
                    Z <- Specs[[2]]
                    if(!is.null(Specs[[2]])) {
                         if(is.matrix(Z)) {
                              if(ncol(Z) != length(Initial.Values))
                                   stop("Z has the wrong number of columns.")
                              if(nrow(Z) != Nc)
                                   stop("Z has the wrong number of rows.")
                              }
                         }
                    }
               if(Specs[["beta"]] == Specs[[3]]) {
                    beta <- Specs[[3]]
                    if(beta <= 1) {
                         cat("\nbeta must be > 1. Changed to 2.\n")
                         beta <- 2}
                    if(length(beta) != 1) {
                         cat("\nLength of beta is wrong. Changed to 1.\n")
                         beta <- beta[1]}}
               else {
                    beta <- 2
                    cat("\nbeta was misspecified and changed to 2.\n")}
               if(Specs[["CPUs"]] == Specs[[4]])
                    CPUs <- max(1, abs(round(Specs[[4]])))
               if(CPUs > 1 & Nc %% 2 != 0)
                    stop("For CPUs > 1, Nc must be even.")
               Packages <- NULL
               if(!is.null(Specs[["Packages"]]))
                    Packages <- Specs[["Packages"]]
               Dyn.libs <- NULL
               if(!is.null(Specs[["Dyn.libs"]]))
                    Dyn.libs <- Specs[["Dyn.libs"]]
               }
          if(Algorithm == "AM") {
               Algorithm <- "Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}}
          if(Algorithm == "AMM") {
               Algorithm <- "Adaptive-Mixture Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               if(!is.null(Specs[["B"]])) {
                    Blocks <- Specs[[2]]
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Blocks)) {
                              Covar[[b]] <- diag(length(Blocks[[b]]))}
                         }
                    }
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[3]]) Periodicity <- Specs[[3]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}
               if(Specs[["w"]] == Specs[[4]]) {
                    w <- Specs[[4]]
                    if(w <= 0 || w >= 1) w <- 0.05}
               else {
                    w <- 0.05
                    cat("\nw was misspecified and changed to 0.05.\n")}}
          if(Algorithm == "AMWG") {
               Algorithm <- "Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 1) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[1]]) Periodicity <- Specs[[1]]
               else {
                    Periodicity <- 50
                    cat("\nPeriodicity was misspecified and changed to 50.\n")}}
          if(Algorithm == "CHARM") {
               Algorithm <- "Componentwise Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs)) {
                    alpha.star <- NA
                    }
               else {
                    if(length(Specs) != 1)
                         stop("The Specs argument is incorrect")
                    if(Specs[["alpha.star"]] == Specs[[1]]) {
                         alpha.star <- abs(Specs[[1]][1])
                         if(alpha.star <= 0 | alpha.star >= 1) {
                              cat("\nalpha.star not in (0,1), set to 0.44.\n")
                              alpha.star <- 0.44}
                         }
                    }
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "DEMC") {
               Algorithm <- "Differential Evolution Markov Chain"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1
               if(Specs[["Nc"]] == Specs[[1]]) Nc <- abs(round(Specs[[1]]))
               if(Nc < 3) Nc <- 3
               if(is.null(Specs[[2]]) | all(Specs[["Z"]] == Specs[[2]])) {
                    Z <- Specs[[2]]
                    if(!is.null(Specs[[2]])) {
                         if(is.matrix(Z)) {
                              if(ncol(Z) != length(Initial.Values))
                                   stop("Z has the wrong number of columns.")
                              if(nrow(Z) != (floor(Iterations/Thinning)+1)) {
                                   Z.temp <- Z[nrow(Z),]
                                   if(nrow(Z) < (floor(Iterations/Thinning)+1)) {
                                        Z <- rbind(Z,
                                             Z[1:(floor(Iterations/Thinning)+1-nrow(Z)),])
                                        }
                                   else if(nrow(Z) > (floor(Iterations/Thinning)+1)) {
                                        Z <- Z[1:(floor(Iterations/Thinning)+1),]
                                        }
                                   Z[1,] <- Z.temp
                                   }
                              Z <- array(Z, dim=c(floor(Iterations/Thinning)+1,
                                   length(Initial.Values), Nc))
                              }
                         if(dim(Z)[1] != floor(Iterations/Thinning)+1)
                              stop("The first dimension of Z is incorrect.")
                         if(dim(Z)[2] != length(Initial.Values))
                              stop("The second dimension of Z is incorrect.")
                         if(dim(Z)[3] != Nc)
                              stop("The third dimension of Z is incorrect.")
                         }
                    }
               if(is.null(Specs[[3]])) {
                    gamma <- 2.381204 / sqrt(2*length(Initial.Values))
                    }
               else if(Specs[["gamma"]] == Specs[[3]]) {
                    gamma <- abs(Specs[[3]])}
               else stop("gamma is misspecified.")
               if(Specs[["w"]] == Specs[[4]]) interval(w <- Specs[[4]], 0, 1)
               }
          if(Algorithm == "DRAM") {
               Algorithm <- "Delayed Rejection Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               DR <- 1
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}}
          if(Algorithm == "DRM") {
               Algorithm <- "Delayed Rejection Metropolis"
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "Experimental") {
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "HARM") {
               Algorithm <- "Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs)) {
                    alpha.star <- NA
                    }
               else {
                    if(length(Specs) != 1)
                         stop("The Specs argument is incorrect")
                    if(Specs[["alpha.star"]] == Specs[[1]]) {
                         alpha.star <- abs(Specs[[1]][1])
                         if(alpha.star <= 0 | alpha.star >= 1) {
                              cat("\nalpha.star not in (0,1), set to 0.234.\n")
                              alpha.star <- 0.234}
                         }
                    }
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "HMC") {
               Algorithm <- "Hamiltonian Monte Carlo"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(all(Specs[["epsilon"]] == Specs[[1]])) {
                    epsilon <- abs(Specs[[1]])
                    if(length(epsilon) != length(Initial.Values))
                         cat("\nLength of epsilon is incorrect.\n")
                         epsilon <- rep(epsilon[1], length(Initial.Values))}
               if(Specs[["L"]] == Specs[[2]]) L <- abs(round(Specs[[2]]))
               if(L < 1) {
                    cat("\nL has been increased to its minimum: 1.\n")
                    L <- 1}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "HMCDA") {
               Algorithm <- "Hamiltonian Monte Carlo with Dual-Averaging"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 5) stop("The Specs argument is incorrect.")
               if(Specs[["A"]] == Specs[[1]])
                    A <- min(round(abs(Specs[[1]])), Iterations)
               if(Specs[["delta"]] == Specs[[2]])
                    delta <- max(min(abs(Specs[[2]]), 1), 1/Iterations)
               if(is.null(Specs[[3]]) |
                    {all(Specs[["epsilon"]] == Specs[[3]])}) {
                    if(is.null(Specs[[3]])) epsilon <- NULL
                    else epsilon <- abs(Specs[[3]][1])}
               if(Specs[["Lmax"]] == Specs[[4]]) 
                    Lmax <- abs(Specs[[4]])
               if(Specs[["lambda"]] == Specs[[5]]) {
                    lambda <- abs(Specs[[5]])
                    if(!is.null(epsilon))
                         if(lambda < epsilon) lambda <- epsilon}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "IM") {
               Algorithm <- "Independence Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 1) stop("The Specs argument is incorrect.")
               if(all(Specs[["mu"]] == Specs[[1]])) {
                    mu <- as.vector(Specs[[1]])
                    if(length(mu) != length(Initial.Values))
                         stop("length(mu) != length(Initial.Values).")}
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "INCA") {
               Algorithm <- "Interchain Adaptation"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}
               }
          if(Algorithm == "MWG") {
               Algorithm <- "Metropolis-within-Gibbs"
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "NUTS") {
               Algorithm <- "No-U-Turn Sampler"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(Specs[["A"]] == Specs[[1]])
                    A <- max(min(round(abs(Specs[[1]])), Iterations),1)
               if(Specs[["delta"]] == Specs[[2]])
                    delta <- max(min(abs(Specs[[2]]), 1), 1/Iterations)
               if(is.null(Specs[[3]]) |
                    {all(Specs[["epsilon"]] == Specs[[3]])}) {
                    if(is.null(Specs[[3]])) epsilon <- NULL
                    else epsilon <- abs(Specs[[3]][1])}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "RAM") {
               Algorithm <- "Robust Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               if(Specs[["alpha.star"]] == Specs[[1]]) {
                    alpha.star <- Specs[[1]]
                    if(alpha.star <= 0 || alpha.star >= 1) {
                         cat("\nalpha.star not in (0,1). Changed to 0.234.\n")
                         alpha.star <- 0.234}
                    if({length(alpha.star) != 1} &
                         {length(alpha.star) != length(Initial.Values)}) {
                         cat("\nLength of alpha.star is wrong. Changed to 1.\n")
                         alpha.star <- alpha.star[1]}}
               else {
                    alpha.star <- 0.234
                    cat("\nalpha.star was misspecified and changed to 0.234.\n")}
               Adaptive <- 2
               DR <- 0
               if(Specs[["Dist"]] == Specs[[2]]) {
                    Dist <- Specs[[2]]
                    if(Dist != "t" & Dist != "N") {
                         cat("\nDist was not t or N, and changed to N.\n")
                         Dist <- "N"}}
               else {
                    Dist <- "N"
                    cat("\nDist was not t or N, and changed to N.\n")
                    }
               if(Specs[["gamma"]] == Specs[[3]]) {
                    gamma <- Specs[[3]]
                    if(gamma <= 0.5 || gamma > 1) {
                         cat("\ngamma not in (0.5,1]. Changed to 0.66.\n")
                         gamma <- 0.66}}
               else {
                    gamma <- 0.66
                    cat("\ngamma was misspecified and changed to 0.66.\n")}
               if(Specs[["Periodicity"]] == Specs[[4]]) Periodicity <- Specs[[4]]
               else {
                    Periodicity <- 10
                    cat("\nPeriodicity was misspecified and changed to 10.\n")}}
          if(Algorithm == "RJ") {
               Algorithm <- "Reversible-Jump"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 5) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1
               if(Specs[["bin.n"]] == Specs[[1]]) bin.n <- round(Specs[[1]])
               else stop("bin.n was misspecified.")
               if(bin.n > length(Initial.Values))
                    bin.n <- length(Initial.Values)
               if(bin.n < 1) bin.n <- 1
               if(Specs[["bin.p"]] == Specs[[2]]) bin.p <- Specs[[2]]
               else stop("bin.p was misspecified.")
               if(bin.p < 0 | bin.p > 1) {
                    bin.p <- interval(bin.p, 0, 1, reflect=FALSE)
                    cat("\nbin.p must be in [0,1]. It's now",
                         round(bin.p,5), "\n")}
               if(all(Specs[["parm.p"]] == Specs[[3]])) parm.p <- Specs[[3]]
               else stop("parm.p was misspecified.")
               if(!is.vector(parm.p)) parm.p <- as.vector(parm.p)
               if(length(parm.p) != length(Initial.Values)) {
                    parm.p <- rep(parm.p[1], length(Initial.Values))
                    cat("\nparm.p now has the correct length, all equal to parm.p[1].\n")}
               if(all(Specs[["selectable"]] == Specs[[4]])) selectable <- Specs[[4]]
               else stop("selectable was misspecified.")
               if(!is.vector(selectable)) selectable <- as.vector(selectable)
               if(length(selectable) != length(Initial.Values)) {
                    selectable <- rep(1, length(Initial.Values))
                    cat("\nselectable now has the correct length, all set to 1.\n")}
               if(all(Specs[["selected"]] == Specs[[5]])) selected <- Specs[[5]]
               else stop("selected was misspecified.")
               if(!is.vector(selected)) selected <- as.vector(selected)
               if(length(selected) != length(Initial.Values)) {
                    selected <- rep(1, length(Initial.Values))
                    cat("\nselected now has the correct length, all set to 1.\n")}}
          if(Algorithm == "RWM") {
               Algorithm <- "Random-Walk Metropolis"
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "SAMWG") {
               Algorithm <- "Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 50
                    cat("\nPeriodicity was misspecified and changed to 50.\n")}}
          if(Algorithm == "Slice") {
               Algorithm <- "Slice Sampler"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(all(Specs[["m"]] == Specs[[1]])) m <- abs(round(Specs[[1]]))
               if(length(m) == 1) m <- rep(m, length(Initial.Values))
               else if(length(m) != length(Initial.Values)) {
                    cat("\nm was misspecified, and is replaced with Inf.\n")
                    m <- rep(Inf, length(Initial.Values))}
               if(any(m < 1)) {
                    cat("\nm was misspecified, and is replaced with 1.\n")
                    m <- ifelse(m < 1, 1, m)}
               if(all(Specs[["w"]] == Specs[[2]])) w <- abs(Specs[[2]])
               if(length(w) == 1) w <- rep(w, length(Initial.Values))
               else if(length(w) != length(Initial.Values)) {
                    cat("\nw was misspecified, and is replaced with 1.\n")
                    w <- rep(1, length(Initial.Values))}
               if(any(w <= 0)) {
                    cat("\nw was misspecified, and is replaced with 1.\n")
                    w <- ifelse(w <= 0, 1, w)}
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "SMWG") {
               Algorithm <- "Sequential Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 1) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)}
          if(Algorithm == "THMC") {
               Algorithm <- "Tempered Hamiltonian Monte Carlo"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(all(Specs[["epsilon"]] == Specs[[1]])) {
                    epsilon <- as.vector(abs(Specs[[1]]))
                    if(length(epsilon) != length(Initial.Values))
                         cat("\nLength of epsilon is incorrect.\n")
                         epsilon <- rep(epsilon[1], length(Initial.Values))}
               if(Specs[["L"]] == Specs[[2]]) L <- abs(round(Specs[[2]]))
               if(L < 2) {
                    cat("\nL has been increased to its minimum: 2.\n")
                    L <- 2}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- 1
               if(Specs[["Temperature"]] == Specs[[3]])
                    Temperature <- Specs[[3]]
                    if(Temperature <= 0) {
                         cat("\nTemperature is incorrect, changed to 1.\n")
                         Temperature <- 1}}
          if(Algorithm == "twalk") {
               Algorithm <- "t-walk"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               if(is.null(Specs[[1]]) | all(Specs[["SIV"]] == Specs[[1]])) {
                    if(is.null(Specs[[1]])) {
                         cat("\nGenerating SIV...\n")
                         if(!is.null(Data$PGF))
                              SIV <- GIV(Model, Data, PGF=TRUE)
                         else SIV <- GIV(Model, Data)}
                    else SIV <- Specs[[1]]
                    if(!identical(length(SIV), length(Initial.Values))) {
                         cat("\nGenerating SIV due to length mismatch.\n")
                         if(!is.null(Data$PGF))
                              SIV <- GIV(Model, Data, PGF=TRUE)
                         else SIV <- GIV(Model, Data)}}
               else {
                    cat("\nSIV was misspecified. Generating values now.\n")
                    if(!is.null(Data$PGF))
                         SIV <- GIV(Model, Data, PGF=TRUE)
                    else SIV <- GIV(Model, Data)}
               Mo2 <- Model(SIV, Data)
               if(!is.finite(Mo2[[1]]))
                    stop("SIV results in a non-finite posterior.")
               if(!is.finite(Mo2[[2]]))
                    stop("SIV results in a non-finite deviance.")
               SIV <- Mo2[[5]]
               rm(Mo2)
               if(Specs[["n1"]] == Specs[[2]]) {
                    n1 <- Specs[[2]]
                    if(n1 < 1) {
                         cat("\nn1 must be at least 1. Changed to 4.\n")
                         n1 <- 4}}
               else {
                    n1 <- 4
                    cat("\nn1 was misspecified and changed to 4.\n")}
               if(Specs[["at"]] == Specs[[3]]) {
                    at <- Specs[[3]]
                    if(at <= 0) {
                         cat("\nat must be positive. Changed to 6.\n")
                         at <- 6}}
               else {
                    at <- 6
                    cat("\nat was misspecified and changed to 6.\n")}
               if(Specs[["aw"]] == Specs[[4]]) {
                    aw <- Specs[[4]]
                    if(aw <= 0) {
                         cat("\naw must be positive. Changed to 1.5.\n")
                         at <- 1.5}}
               else {
                    aw <- 1.5
                    cat("\naw was misspecified and changed to 1.5.\n")}
               Periodicity <- 1}
          if(Algorithm == "USAMWG") {
               Algorithm <- "Updating Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 50
                    cat("\nPeriodicity was misspecified and changed to 50.\n")}
               if({class(Specs[["Fit"]]) == "demonoid"} &
                    {class(Specs[[3]]) == "demonoid"}) Fit <- Specs[[3]]
               if(Specs[["Begin"]] == Specs[[4]]) Begin <- Specs[[4]]}
          if(Algorithm == "USMWG") {
               Algorithm <- "Updating Sequential Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               Periodicity <- Iterations + 1
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)
               if({class(Specs[["Fit"]]) == "demonoid"} &
                    {class(Specs[[2]]) == "demonoid"}) Fit <- Specs[[2]]
               if(Specs[["Begin"]] == Specs[[3]]) Begin <- Specs[[3]]}
               }
     else{
          cat("Unknown algorithm has been changed to Random-Walk Metropolis.\n")
          Algorithm <- "Random-Walk Metropolis"
          Adaptive <- Iterations + 1
          DR <- 0
          Periodicity <- Iterations + 1
          }
     #     
     Adaptive <- round(Adaptive)
     if({Adaptive < 1} || {Adaptive > Iterations}) 
          Adaptive <- Iterations + 1
     Periodicity <- round(Periodicity)
     if({Periodicity < 1} || {Periodicity > Iterations}) 
          Periodicity <- Iterations + 1
     Mo0 <- Model(Initial.Values, Data)
     if(!is.list(Mo0)) stop("Model must return a list.")
     if(length(Mo0) != 5) stop("Model must return five components.")
     if(any(names(Mo0) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.")
     if(length(Mo0[["LP"]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(Mo0[["Monitor"]]), length(Data$mon.names)))
          stop("Length of mon.names differs from length of monitors.")
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out, 1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of apply functions\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are 'vectorized'.\n")}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are 'vectorized'.\n")}
#############################################################################

     #########################  Initial Settings  #########################
     Acceptance <- 0
     if(!is.finite(Mo0[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          if(!is.null(Data$PGF))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data)
          Mo0 <- Model(Initial.Values, Data)
          }
     if(is.infinite(Mo0[["LP"]])) stop("The posterior is infinite!")
     if(is.nan(Mo0[["LP"]])) stop("The posterior is not a number!")
     if(is.na(Mo0[["Dev"]])) stop("The deviance is a missing value!")
     if(is.infinite(Mo0[["Dev"]])) stop("The deviance is infinite!")
     if(is.nan(Mo0[["Dev"]])) stop("The deviance is not a number!")
     if(any(is.na(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have a missing value!")
     if(any(is.infinite(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have an infinite value!")
     if(any(is.nan(Mo0[["Monitor"]])))
          stop("Monitored variable(s) include a value that is not a number!")
     if(Algorithm == "t-walk") {
          Mo0 <- Model(Initial.Values, Data)
          if(any(Mo0[["parm"]] == SIV))
              stop("Initial.Values and SIV not unique after model update.")}
     ######################  Laplace Approximation  #######################
     ### Sample Size of Data
     if(!is.null(Data$n)) if(length(Data$n) == 1) N <- Data$n
     if(!is.null(Data$N)) if(length(Data$N) == 1) N <- Data$N
     if(!is.null(Data$y)) N <- nrow(matrix(Data$y))
     if(!is.null(Data$Y)) N <- nrow(matrix(Data$Y))
     if(is.null(N)) stop("Sample size of Data not found in n, N, y, or Y.")
     if({all(Initial.Values == 0)} & {N >= 5*length(Initial.Values)}) {
          cat("\nLaplace Approximation will be used on initial values.\n")
          Fit.LA <- LaplaceApproximation(Model, Initial.Values, Data,
               Method="HAR", sir=FALSE)
          Covar <- 2.381204 * 2.381204 / length(Initial.Values) *
               Fit.LA$Covar
          Initial.Values <- Fit.LA$Summary1[1:length(Initial.Values),1]
          cat("The covariance matrix from Laplace Approximation has been scaled\n")
          cat("for Laplace's Demon, and the posterior modes are now the initial\n")
          cat("values for Laplace's Demon.\n\n")}
     #########################  Prepare for MCMC  #########################
     Mo0 <- Model(Initial.Values, Data)
     Dev <- matrix(Mo0[["Dev"]], floor(Iterations/Thinning)+1, 1)
     Mon <- matrix(Mo0[["Monitor"]], floor(Iterations/Thinning)+1,
          length(Mo0[["Monitor"]]), byrow=TRUE)
     LIV <- length(Initial.Values)
     thinned <- matrix(0, floor(Iterations/Thinning)+1,
          length(Initial.Values))
     thinned[1,] <- prop <- Initial.Values
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(Algorithm %in% c("Adaptive Metropolis",
          "Adaptive-Mixture Metropolis",
          "Delayed Rejection Adaptive Metropolis",
          "Delayed Rejection Metropolis", "Interchain Adaptation",
          "Random-Walk Metropolis")) {
          ### Algorithms that require both VarCov and tuning
  #        
          if(is.list(Covar) & Algorithm != "Adaptive-Mixture Metropolis") {
               Covar <- NULL}
          else if(is.matrix(Covar) & !is.list(Covar)) {
               tuning <- sqrt(diag(Covar)); VarCov <- Covar}
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != LIV)
                    tuning <- rep(ScaleF, LIV)
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning
               }
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning}
          else if(is.list(Covar)) {
               tuning <- Covar
               for (i in 1:length(tuning)) {
                    tuning[[i]] <- sqrt(diag(tuning[[i]]))}
               VarCov <- Covar}
          if(is.matrix(VarCov) & !is.list(VarCov)) {
               DiagCovar <- matrix(diag(VarCov), 1, LIV)}
          else if(is.list(VarCov)) {
               DiagCovar <- matrix(1, 1, LIV)
               for (b in 1:length(Blocks)) {
                    DiagCovar[Blocks[[b]]] <- diag(VarCov[[b]])}}
          }
     else if(Algorithm %in% c("Independence Metropolis",
          "Robust Adaptive Metropolis")) {
          ### Algorithms that require VarCov, but not tuning
          if(is.list(Covar)) {
               Covar <- NULL}
          else if(is.matrix(Covar) & !is.list(Covar)) {
               VarCov <- Covar}
          else if(is.vector(Covar) & !is.list(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- abs(as.vector(Covar))
               }
          else if(is.null(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- rep(ScaleF, LIV)}
          else if(is.list(Covar)) {VarCov <- Covar}
          if(is.matrix(VarCov) & !is.list(VarCov)) {
               DiagCovar <- matrix(diag(VarCov), 1, LIV)}
          else if(is.list(VarCov)) {
               DiagCovar <- matrix(1, 1, LIV)
               for (b in 1:length(Blocks)) {
                    DiagCovar[Blocks[[b]]] <- diag(VarCov[[b]])}}
          }
     else if(Algorithm %in% c("Adaptive Metropolis-within-Gibbs",
          "Metropolis-within-Gibbs",
          "Sequential Adaptive Metropolis-within-Gibbs",
          "Sequential Metropolis-within-Gibbs",
          "Updating Sequential Adaptive Metropolis-within-Gibbs",
          "Updating Sequential Metropolis-within-Gibbs")) {
          ### Algorithms that do not require VarCov, but require tuning
          if(is.list(Covar)) {
               Covar <- NULL}
          else if(is.matrix(Covar) & !is.list(Covar)) {
               tuning <- sqrt(diag(Covar))}
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != length(Initial.Values))
                    tuning <- rep(ScaleF, LIV)}
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)}
          VarCov <- NULL
          DiagCovar <- matrix(tuning, 1, LIV)
          }
     else {
          ### Algorithms that do not require VarCov or tuning
          VarCov <- NULL
          DiagCovar <- matrix(1, 1, LIV)
          }
     rm(Covar)
     ############################  Begin MCMC  ############################
     cat("Algorithm:", Algorithm, "\n")
     cat("\nLaplace's Demon is beginning to update...\n")
     if(Algorithm == "Adaptive Hamiltonian Monte Carlo") {
          mcmc.out <- AHMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, epsilon, L)}
     else if(Algorithm == "Affine-Invariant Ensemble Sampler") {
          mcmc.out <- AIES(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, Nc, Z, beta, CPUs,
               Packages, Dyn.libs)}
     else if(Algorithm == "Adaptive Metropolis") {
          mcmc.out <- AM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & !is.list(VarCov)) {
          mcmc.out <- AMM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, w)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & is.list(VarCov)) {
          mcmc.out <- AMM.B(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, Blocks, w)}
     else if(Algorithm == "Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- AMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning)}
     else if(Algorithm == "Componentwise Hit-And-Run Metropolis") {
          mcmc.out <- CHARM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, alpha.star)}
     else if(Algorithm == "Delayed Rejection Adaptive Metropolis") {
          mcmc.out <- DRAM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov)}
     else if(Algorithm == "Delayed Rejection Metropolis") {
          mcmc.out <- DRM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov)}
     else if(Algorithm == "Differential Evolution Markov Chain") {
          mcmc.out <- DEMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, Nc, Z, gamma, w)}
     else if(Algorithm == "Experimental") {
#          mcmc.out <- Experimental(Model, Data, Adaptive, DR, Iterations,
#               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
#               LIV, Mon, Mo0, ScaleF, thinned)}
          stop("Experimental function not found.")}
     else if(Algorithm == "Hamiltonian Monte Carlo") {
          mcmc.out <- HMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, epsilon, L)}
     else if(Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          mcmc.out <- HMCDA(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, A, delta, epsilon, Lmax,
               lambda)}
     else if(Algorithm == "Hit-And-Run Metropolis") {
          mcmc.out <- HARM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, alpha.star)}
     else if(Algorithm == "Independence Metropolis") {
          mcmc.out <- IM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, VarCov, mu)}
     else if(Algorithm == "Interchain Adaptation") {
          mcmc.out <- INCA(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov)}
     else if(Algorithm == "Metropolis-within-Gibbs") {
          mcmc.out <- MWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning)}
     else if(Algorithm == "No-U-Turn Sampler") {
          mcmc.out <- NUTS(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, A, delta, epsilon)}
     else if(Algorithm == "Robust Adaptive Metropolis") {
          mcmc.out <- RAM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, VarCov, alpha.star, Dist,
               gamma)}
     else if(Algorithm == "Random-Walk Metropolis") {
          mcmc.out <- RWM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov)}
     else if(Algorithm == "Reversible-Jump") {
          mcmc.out <- RJ(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, bin.n, bin.p, parm.p,
               selectable, selected)}
     else if(Algorithm == "Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- SAMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning,
               parm.names=Data$parm.names, Dyn)}
     else if(Algorithm == "Sequential Metropolis-within-Gibbs") {
          mcmc.out <- SMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning,
               parm.names=Data$parm.names, Dyn)}
     else if(Algorithm == "Slice Sampler") {
          mcmc.out <- Slice(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, m, w)}
     else if(Algorithm == "Tempered Hamiltonian Monte Carlo") {
          mcmc.out <- THMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, epsilon, L, Temperature)}
     else if(Algorithm == "t-walk") {
          mcmc.out <- twalk(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, SIV=SIV, n1=n1, at=at,
               aw=aw)}
     else if(Algorithm == "Updating Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- USAMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning,
               parm.names=Data$parm.names, Dyn, Fit, Begin)}
     else if(Algorithm == "Updating Sequential Metropolis-within-Gibbs") {
          mcmc.out <- USMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               LIV, Mon, Mo0, ScaleF, thinned, tuning,
               parm.names=Data$parm.names, Dyn, Fit, Begin)}
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


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#

 
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################################################################################
AHMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, epsilon, L)
     {
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     DiagCovar[1,] <- epsilon
     gr0 <- partial(Model, post[1,], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- post[iter,]
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
               post[iter,] <- Mo1[["parm"]]
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Adaptation
          if({iter > 10} & {iter %% Periodicity == 0}) {
               acceptances <- apply(post[(iter-9):iter,], 2, function(x)
                    {length(unique(x))})
               eps.num <- which(acceptances <= 1)
               epsilon[eps.num] <- epsilon[eps.num] * 0.8
               eps.num <- which(acceptances > 7)
               epsilon[eps.num] <- epsilon[eps.num] * 1.2
               DiagCovar <- rbind(DiagCovar, epsilon)}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
AIES <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, Nc, Z, beta, CPUs, Packages, Dyn.libs)
     {
     Mo0 <- list(Mo0=Mo0)
     if(is.null(Z)) {
          Z <- matrix(Mo0[[1]][["parm"]], Nc, LIV, byrow=TRUE)
          for (i in 2:Nc) {
               if(!is.null(Data$PGF)) {
                    Z[i,] <- GIV(Model, Data, PGF=TRUE)
                    }
               else Z[i,] <- GIV(Model, Data)
               }
          }
     for (i in 2:Nc) Mo0[[i]] <- Model(Z[i,], Data)
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[[1]][["parm"]]
                    Dev[t.iter] <- Mo0[[1]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
               for (i in 1:Nc) {
                    ### Propose new parameter values with stretch move
                    z <- 1 / sqrt(runif(1, 1 / beta, beta))
                    s <- sample(c(1:Nc)[-i], 1)
                    prop <- Mo0[[s]][["parm"]] +
                         z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
                    if(i == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Multivariate\n")
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0[[i]]
                    if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0[[i]]
                    if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0[[i]]
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[["parm"]]
                                   Dev[t.iter] <- Mo1[["Dev"]]
                                   Mon[t.iter,] <- Mo1[["Monitor"]]}
                              }
                         }
                    }
               }
          }
     else {
          library(parallel, quietly=TRUE)
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores,"\n")
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n")
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...")
          cat("\n##################################################\n")
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n")
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
          ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
               envir=environment())
          model.wrapper <- function(x, ...)
               {
               if(!is.null(Packages)) {
                    sapply(Packages,
                         function(x) library(x, character.only=TRUE,
                              quietly=TRUE))}
               if(!is.null(Dyn.libs)) {
                    sapply(Dyn.libs,
                         function(x) dyn.load(paste(wd, x, sep = "/")))
                    on.exit(sapply(Dyn.libs,
                         function(x) dyn.unload(paste(wd, x, sep = "/"))))}
               Model(prop[x,], Data)
               }
          prop <- Z
          batch1 <- 1:(Nc/2)
          batch2 <- batch1 + (Nc/2)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[[1]][["parm"]]
                    Dev[t.iter] <- Mo0[[1]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
               for (i in 1:Nc) {
                    ### Propose new parameter values with stretch move
                    z <- 1 / sqrt(runif(1, 1 / beta, beta))
                    if(i <= (Nc/2)) s <- sample(batch2, 1)
                    else s <- sample(batch1, 1)
                    prop[i,] <- Mo0[[s]][["parm"]] +
                         z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
                    if(i == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Multivariate\n")}
               ### Log-Posterior of the proposed state
               Mo1 <- clusterApply(cl, 1:Nc, model.wrapper,
                    Model, Data, prop)
               for (i in 1:Nc) {
                    if(!is.finite(Mo1[[i]][["LP"]])) Mo1[[i]] <- Mo0[[i]]
                    if(!is.finite(Mo1[[i]][["Dev"]])) Mo1[[i]] <- Mo0[[i]]
                    if(any(!is.finite(Mo1[[i]][["Monitor"]])))
                         Mo1[[i]] <- Mo0[[i]]
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[[i]][["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1[[i]]
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[[i]][["parm"]]
                                   Dev[t.iter] <- Mo1[[i]][["Dev"]]
                                   Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}
                              }
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
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
AMM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, VarCov, w)
     {
     obs.sum <- matrix(0, LIV, 1)
     obs.scatter <- matrix(0, LIV, LIV)
     prop.R <- NULL
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values from a mixture
          if(is.null(prop.R) || runif(1) < w) {
               prop <- rnorm(LIV, Mo0[["parm"]], tuning)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Non-Adaptive Component\n")}
          else {
               prop <- t(prop.R) %*% rnorm(LIV) + Mo0[["parm"]]
               if(iter %% Status == 0) 
                    cat(",   Proposal: Adaptive Component\n")}
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
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               VarCov <- obs.scatter/iter - tcrossprod(obs.sum/iter)
               diag(VarCov) <- diag(VarCov) + 1e-05
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               prop.R <- try(ScaleF * chol(VarCov), silent=TRUE)
               if(!is.matrix(prop.R)) prop.R <- NULL}
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
AMM.B <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, VarCov, Blocks, w)
     {
     B <- length(Blocks)
     obs.scatter <- obs.sum <- prop.R <- list(NULL)
     for (b in 1:B) {
          obs.sum[[b]] <- matrix(0, length(Blocks[[b]]), 1)
          obs.scatter[[b]] <- matrix(0, length(Blocks[[b]]), length(Blocks[[b]]))
          prop.R[[b]] <- list(NA)
          }
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Proceed by Block
          for (b in 1:B) {
               ### Propose new parameter values from a mixture
               prop <- Mo0[["parm"]]
               if(is.na(prop.R[[b]]) || runif(1) < w) {
                    prop[Blocks[[b]]] <- rnorm(length(Blocks[[b]]),
                         Mo0[["parm"]][Blocks[[b]]], tuning)
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Non-Adaptive Component\n")}
               else {
                    prop[Blocks[[b]]] <- t(prop.R[[b]]) %*%
                         rnorm(length(Blocks[[b]])) +
                         Mo0[["parm"]][Blocks[[b]]]
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Adaptive Component\n")}
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
                    Acceptance <- Acceptance + length(Blocks[[b]]) / LIV
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               ### Update Sample and Scatter Sum
               obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Blocks[[b]]]
               obs.scatter[[b]] <- obs.scatter[[b]] +
                    tcrossprod(Mo0[["parm"]][Blocks[[b]]])
               ### Adapt the Proposal Variance
               if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
                    VarCov[[b]] <- obs.scatter[[b]]/iter -
                         tcrossprod(obs.sum[[b]]/iter)
                    diag(VarCov[[b]]) <- diag(VarCov[[b]]) + 1e-05
                    if(b == 1) DiagCovar <- rbind(DiagCovar, rep(0,LIV))
                    DiagCovar[nrow(DiagCovar),Blocks[[b]]] <- diag(VarCov[[b]])
                    prop.R[[b]] <- try(ScaleF * chol(VarCov[[b]]), silent=TRUE)
                    if(!is.matrix(prop.R[[b]])) prop.R[[b]] <- NA}
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
AMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
CHARM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, alpha.star)
     {
     if(is.na(alpha.star)) {
          Acceptance <- matrix(0, 1, LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) {cat("Iteration: ", iter,
                    ",   Proposal: Componentwise\n", sep="")}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
                    if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
                    if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + u}
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          ### Output
          out <- list(Acceptance=mean(as.vector(Acceptance)),
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     else {
          tau <- rep(1, LIV)
          Acceptance <- matrix(0, 1, LIV)
          DiagCovar <- matrix(tau, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) {cat("Iteration: ", iter,
                    ",   Proposal: Componentwise\n", sep="")}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + tau[j]*lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- Model(prop, Data)
                    if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
                    if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
                    if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         tau[j] <- tau[j] + (tau[j] / (alpha.star *
                         (1 - alpha.star))) * (1 - alpha.star) / iter
                         }
                    else {
                         tau[j] <- abs(tau[j] - (tau[j] / (alpha.star *
                         (1 - alpha.star))) * alpha.star / iter)
                         }
                    Acceptance[j] <- Acceptance[j] + u}
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]
                    DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=mean(as.vector(Acceptance)),
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     }
DEMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, Nc, Z, gamma, w)
     {
     const <- 2.381204 / sqrt(2)
     Mo0 <- list(Mo0=Mo0)
     if(is.null(Z)) {
          cat("\nGenerating Z...\n")
          Z <- array(0, dim=c(floor(Iterations/Thinning)+1, LIV, Nc))
          for (t in 1:dim(Z)[1]) {
               for (i in 1:Nc) {
                    if(t == 1 & i == 1) {
                         Z[t,,i] <- Mo0[[1]][["parm"]]
                         }
                    else {
                         if(!is.null(Data$PGF)) {
                              Z[t,,i] <- GIV(Model, Data, PGF=TRUE)}
                         else Z[t,,i] <- GIV(Model, Data)
                         }
                    }
               }
          }
     else Z[1,,1] <- Mo0[[1]][["parm"]]
     for (i in 2:Nc) Mo0[[i]] <- Model(Z[1,,i], Data)
     for (iter in 1:Iterations) {
          ### Thinned Iteration
          t.iter <- floor(iter / Thinning) + 1
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               Z[t.iter,,] <- Z[t.iter-1,,]
               thinned[t.iter,] <- Mo0[[1]][["parm"]]
               Dev[t.iter] <- Mo0[[1]][["Dev"]]
               Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
          omega <- runif(1)
          for (i in 1:Nc) {
               r <- sample(dim(Z)[1], 2)
               s <- sample(c(1:Nc)[-i], 2)
               if(omega > w) {
                    ### Parallel Direction Move
                    prop <- Mo0[[i]][["parm"]] +
                         gamma*(Z[r[1],,s[1]] - Z[r[2],,s[2]]) +
                         runif(LIV, -0.001, 0.001)^LIV
                         }
               else {
                    ### Snooker Move
                    si <- sample(c(1:Nc)[-i], 1)
                    prop <- Mo0[[i]][["parm"]] + const*
                         ({Mo0[[si]][["parm"]] - Z[r[1],,s[1]]} -
                          {Mo0[[si]][["parm"]] - Z[r[2],,s[2]]})}
               if(i == 1 & iter %% Status == 0) 
                    cat(",   Proposal: Multivariate\n")
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0[[i]]
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0[[i]]
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0[[i]]
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[[i]][["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0[[i]] <- Mo1
                    Z[t.iter,,i] <- Mo1[["parm"]]
                    if(i == 1) {
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=thinned,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
DRAM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
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
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov), silent=TRUE)
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
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVN.rand <- rnorm(LIV, 0, 1)
               MVNz <- try(matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- Model(prop, Data)
               if(!is.finite(Mo2[["LP"]])) Mo2 <- Mo0
               if(!is.finite(Mo2[["Dev"]])) Mo2 <- Mo0
               if(any(!is.finite(Mo2[["Monitor"]]))) Mo2 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[["LP"]] - Mo2[["LP"]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                    {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] - Mo0[["LP"]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo2
                    post[iter,] <- Mo2[["parm"]]
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
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
DRM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
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
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVN.rand <- rnorm(LIV, 0, 1)
               MVNz <- try(matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
               else {
                    prop <- Mo0[["parm"]]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- Model(prop, Data)
               if(!is.finite(Mo2[["LP"]])) Mo2 <- Mo0
               if(!is.finite(Mo2[["Dev"]])) Mo2 <- Mo0
               if(any(!is.finite(Mo2[["Monitor"]]))) Mo2 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[["LP"]] - Mo2[["LP"]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                    {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] - Mo0[["LP"]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo2
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
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
HARM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, alpha.star)
     {
     if(is.na(alpha.star)) {
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
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1) * d
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate\n")
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
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     else {
          tau <- 1
          DiagCovar <- matrix(tau, nrow(thinned), LIV)
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
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1,0,tau) * d
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate\n")
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
                    tau <- tau + (tau / (alpha.star *
                         (1 - alpha.star))) * (1 - alpha.star) / iter
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]
                         DiagCovar[t.iter,] <- tau}
                    }
               else {
                    tau <- abs(tau - (tau / (alpha.star *
                         (1 - alpha.star))) * alpha.star / iter)
                    if(iter %% Thinning == 0) DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
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
HMCDA <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, A, delta, epsilon, Lmax, lambda)
     {
     leapfrog <- function(theta, r, grad, epsilon, Model, Data)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- Model(thetaprime, Data)
          thetaprime <- Mo1[["parm"]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime,
               rprime=rprime,
               gradprime=gradprime,
               Mo1=Mo1)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data)
          {
          cat("\nFinding a reasonable initial value for epsilon...")
          epsilon <- 1
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
          if(!is.finite(leap$Mo1[["LP"]]))
               stop("LP is not finite in find.reasonable.epsilon().")
          acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
               if(!is.finite(leap$Mo1[["LP"]]))
                    stop("LP is not finite in find.reasonable.epsilon().")
               acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="")
          return(epsilon)
          }
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(Mo0[["parm"]], gr0, Mo0, Model,
               Data)
     DiagCovar[1,] <- epsilon
     L <- max(1, round(lambda / epsilon))
     L <- min(L, Lmax)
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Begin HMCDA
     for (iter in 1:Iterations) {
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- runif(LIV)
          joint <- Mo0[["LP"]] - 0.5 * as.vector(momentum0 %*% momentum0)
          L <- max(1, round(lambda / epsilon))
          L <- min(L, Lmax)
          gr1 <- gr0
          Mo0.1 <- Mo0
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate,   L: ", L, "\n", sep="")
          ### Leapfrog Function
          for (l in 1:L) {
               momentum1 <- momentum1 + 0.5 * epsilon * gr1
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
               momentum1 <- momentum1 + epsilon * gr1}
          ### Accept/Reject
          alpha <- min(1,
               exp(prop - 0.5 * as.vector(momentum1 %*% momentum1) - joint))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Adaptation
          if(iter > 1) {
               eta <- 1 / (iter - 1 + t0)
               Hbar <- (1 - eta) * Hbar + eta * (delta - alpha)
               if(iter <= A) {
                    epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
                    eta <- (iter-1)^-kappa
                    epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                         eta * log(epsilon))
                    DiagCovar <- rbind(DiagCovar, epsilon)}
               else epsilon <- epsilonbar}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
IM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, VarCov, mu)
     {
     VarCov2 <- as.positive.definite(as.symmetric.matrix(VarCov * 1.1))
     Omega <- as.inverse(VarCov2)
     d <- eigen(VarCov2, symmetric=TRUE)$values
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
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov2),
               silent=TRUE)
          if(!inherits(MVNz, "try-error")) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- as.vector(t(mu + t(MVNz)))}
          else {prop <- as.vector(Mo0[["parm"]])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
          if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
          ### Importance Densities (dmvn)
          ss <- prop - mu
          z <- rowSums({ss %*% Omega} * ss)
          d1 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
          ss <- Mo0[["parm"]] - mu
          z <- rowSums({ss %*% Omega} * ss)
          d0 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + d1 - d0
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
          VarCov=cov(thinned))
     return(out)
     }
INCA <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, VarCov)
     {
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     con <- get("con")
     Chains <- get("Chains")
     ### Store all posteriors
     INCA_first <- TRUE
     INCA_iter <- 0
     INCA_post <- matrix(NA, nrow=Iterations*Chains, ncol=LIV)
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
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               select_post <- post[(iter-Periodicity+1):iter,]
               ### Ask for last posteriors to hpc_server
               tmp <- unserialize(con)
               ### Send new posteriors matrix to hpc_server      
               serialize(select_post, con)
               for (i in 1:length(select_post[, 1])) {
                    INCA_post[(INCA_iter + i), ] <- select_post[i, ]}
               INCA_iter <- INCA_iter + length(select_post[, 1])
               if(is.matrix(tmp) && INCA_first == FALSE) {
                    for (i in 1:length(tmp[, 1])) {
                         INCA_post[(INCA_iter + i), ] <- tmp[i, ]}
                    INCA_iter <- INCA_iter + length(tmp[, 1])}
               else if(INCA_first == TRUE) INCA_first <- FALSE
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(INCA_post[1:INCA_iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF *
                         {var(INCA_post[1:INCA_iter,j])} + ScaleF * 1.0E-5)}}
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
MWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
NUTS <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, A, delta, epsilon)
     {
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     leapfrog <- function(theta, r, grad, epsilon, Model, Data)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- Model(thetaprime, Data)
          thetaprime <- Mo1[["parm"]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime,
               rprime=rprime,
               gradprime=gradprime,
               Mo1=Mo1)
          return(out)
          }
     stop.criterion <- function(thetaminus, thetaplus, rminus, rplus)
          {
          thetavec <- thetaplus - thetaminus
          criterion <- (thetavec %*% rminus >= 0) &&
               (thetavec %*% rplus >= 0)
          return(criterion)
          }
     build.tree <- function(theta, r, grad, logu, v, j, epsilon, joint0)
          {
          if(j == 0) {
               ### Base case: Take a single leapfrog step in direction v
               leap <- leapfrog(theta, r, grad, v*epsilon, Model, Data)
               rprime <- leap$rprime
               thetaprime <- leap$thetaprime
               Mo1 <- leap$Mo1
               gradprime <- leap$gradprime
               joint <- Mo1[["LP"]] - 0.5 * as.vector(rprime %*% rprime)
               ### Is the new point in the slice?
               nprime <- logu < joint
               ### Is the simulation wildly inaccurate?
               sprime <- logu - 1000 < joint
               # Set the return values---minus=plus for all things here,
               # since the "tree" is of depth 0.
               thetaminus <- thetaprime
               thetaplus <- thetaprime
               rminus <- rprime
               rplus <- rprime
               gradminus <- gradprime
               gradplus <- gradprime
               ### Compute the acceptance probability
               alphaprime <- min(1, exp(Mo1[["LP"]] - 0.5 *
                    as.vector(rprime %*% rprime) - joint0))
               nalphaprime <- 1}
          else {
               # Recursion: Implicitly build the height j-1 left and
               # right subtrees
               tree <- build.tree(theta, r, grad, logu, v, j-1, epsilon,
                    joint0)
               thetaminus <- tree$thetaminus
               rminus <- tree$rminus
               gradminus <- tree$gradminus
               thetaplus <- tree$thetaplus
               rplus <- tree$rplus
               gradplus <- tree$gradplus
               thetaprime <- tree$thetaprime
               gradprime <- tree$gradprime
               Mo1 <- tree$Mo1
               nprime <- tree$nprime
               sprime <- tree$sprime
               alphaprime <- tree$alphaprime
               nalphaprime <- tree$nalphaprime
               ### If the first subtree stopping criterion is met, then stop
               if(sprime == 1) {
                    if(v == -1) {
                         tree <- build.tree(thetaminus, rminus, gradminus,
                              logu, v, j-1, epsilon, joint0)
                         thetaminus <- tree$thetaminus
                         rminus <- tree$rminus
                         gradminus <- tree$gradminus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    else {
                         tree <- build.tree(thetaplus, rplus, gradplus,
                              logu, v, j-1, epsilon, joint0)
                         thetaplus <- tree$thetaplus
                         rplus <- tree$rplus
                         gradplus <- tree$gradplus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    ### Choose a subtree to propagate a sample up from
                    temp <- nprime2 / (nprime + nprime2)
                    if(!is.finite(temp)) temp <- 0
                    if(runif(1) < temp) {
                         thetaprime <- thetaprime2
                         gradprime <- gradprime2
                         Mo1 <- Mo12}
                    ### Update the number of valid points
                    nprime <- nprime + nprime2
                    ### Update the stopping criterion
                    sprime <- sprime && sprime2 &&
                         stop.criterion(thetaminus, thetaplus, rminus,
                              rplus)
                    ### Update acceptance probability statistics
                    alphaprime <- alphaprime + alphaprime2
                    nalphaprime <- nalphaprime + nalphaprime2}}
          out <- list(thetaminus=thetaminus,
               rminus=rminus,
               gradminus=gradminus,
               thetaplus=thetaplus,
               rplus=rplus,
               gradplus=gradplus,
               thetaprime=thetaprime,
               gradprime=gradprime,
               Mo1=Mo1,
               nprime=nprime,
               sprime=sprime,
               alphaprime=alphaprime,
               nalphaprime=nalphaprime)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data)
          {
          cat("\nFinding a reasonable initial value for epsilon...")
          epsilon <- 1
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
          if(!is.finite(leap$Mo1[["LP"]]))
               stop("LP is not finite in find.reasonable.epsilon().")
          acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
               if(!is.finite(leap$Mo1[["LP"]]))
                    stop("LP is not finite in find.reasonable.epsilon().")
               acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="")
          return(epsilon)
          }
     Count <- 0
     evals <- 0
     grad <- partial(Model, post[1,], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(post[1,], grad, Mo0, Model,
               Data)
     DiagCovar[1,] <- epsilon
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Reset Dev, Mon, and thinned
     if(A < Iterations) {
          Dev <- matrix(Dev[1:(floor((Iterations-A)/Thinning)+1),])
          Mon <- matrix(Mo0[["Monitor"]], floor((Iterations-A)/Thinning)+1,
               length(Mo0[["Monitor"]]), byrow=TRUE)
          thinned <- matrix(0, floor((Iterations-A)/Thinning)+1, LIV)}
     ### Begin NUTS
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter > A) {
               if((iter-A) %% Thinning == 0) {
                    thinned[((iter-A)/Thinning+1),] <- post[iter,]
                    Dev[((iter-A)/Thinning+1)] <- Mo0[["Dev"]]
                    Mon[((iter-A)/Thinning+1),] <- Mo0[["Monitor"]]}}
          else if(A >= Iterations) {
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- post[iter,]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}}
          prop <- post[iter,]
          r0 <- runif(LIV) ### r0 is momenta
          ### Joint log-probability of theta and momenta r
          joint <- Mo0[["LP"]] - 0.5 * as.vector(r0 %*% r0)
          ### Resample u ~ U([0, exp(joint)])
          logu <- joint - rexp(1)
          ### Initialize Tree
          thetaminus <- prop
          thetaplus <- prop
          rminus <- r0
          rplus <- r0
          gradminus <- grad
          gradplus <- grad
          j <- 0 ### Initial height j=0
          n <- 1 ### Initially, the only valid point is the initial point
          s <- 1 ### Loop until s == 0
          while (s == 1) {
               ### Choose a direction: -1=backwards, 1=forwards.
               v <- 2*(runif(1) < 0.5) - 1
               ### Double the size of the tree.
               if(v == -1) {
                    tree <- build.tree(thetaminus, rminus, gradminus,
                         logu, v, j, epsilon, joint)
                    thetaminus <- tree$thetaminus
                    rminus <- tree$rminus
                    gradminus <- tree$gradminus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               else {
                    tree <- build.tree(thetaplus, rplus, gradplus, logu,
                         v, j, epsilon, joint)
                    thetaplus <- tree$thetaplus
                    rplus <- tree$rplus
                    gradplus <- tree$gradplus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               ### Accept/Reject
               Count <- Count + 1
               if((sprime == 1) && (runif(1) < nprime/n)) {
                    post[iter,] <- thetaprime
                    Mo0 <- Mo1
                    grad <- gradprime
                    Acceptance <- Acceptance + 1
                    if(iter > A) {
                         if((iter-A) %% Thinning == 0) {
                              thinned[((iter-A)/Thinning+1),] <- Mo1[["parm"]]
                              Dev[((iter-A)/Thinning+1)] <- Mo1[["Dev"]]
                              Mon[((iter-A)/Thinning+1),] <- Mo1[["Monitor"]]}}
                    else if(A >= Iterations) {
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}}}
               ### Update number of observed valid points
               n <- n + nprime
               ### Decide if it is time to stop
               s <- sprime &&
                    stop.criterion(thetaminus, thetaplus, rminus, rplus)
               ### Increment depth
               j <- j + 1}
          ### Adaptation of epsilon
          eta <- 1 / (iter - 1 + t0)
          Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
          if(iter <= A) {
               epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
               eta <- (iter-1)^-kappa
               epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                    eta * log(epsilon))
               DiagCovar <- rbind(DiagCovar, epsilon)}
          else epsilon <- epsilonbar
          }
     Acceptance <- round(Acceptance / Count * Iterations)
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
RAM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, VarCov, alpha.star, Dist, gamma)
     {
     if(!is.symmetric.matrix(VarCov)) {
          cat("\nAsymmetric VarCov, correcting now...\n")
          VarCov <- as.symmetric.matrix(VarCov)}
     if(!is.positive.definite(VarCov)) {
          cat("\nNon-Positive-Definite VarCov, correcting now...\n")
          VarCov <- as.positive.definite(VarCov)}
     Iden.Mat <- diag(LIV)
     S.z <- try(t(chol(VarCov)), silent=TRUE)
     if(!inherits(S.z, "try-error")) S <- S.z
     else S <- Iden.Mat
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose New Parameter Values
          if(Dist == "t") U <- qt(runif(LIV), df=5, lower.tail=TRUE)
          else U <- rnorm(LIV)
          prop <- Mo0[["parm"]] + S %*% U
          if(iter %% Status == 0)
               cat(",   Proposal: Multivariate\n")
          ### Log-Posterior
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
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
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Adaptation
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               eta <- min(1, LIV*iter^(-gamma))
               VarCov.test <- S %*% (Iden.Mat +
                    eta*(min(1, exp(log.alpha)) - alpha.star) *
                    U %*% t(U) / sqrt(sum(U^2))) %*% t(S)
               if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
                    !is.matrix(VarCov.test)) {VarCov.test <- VarCov}
               if(!is.symmetric.matrix(VarCov.test))
                    VarCov.test <- as.symmetric.matrix(VarCov.test)
               if(is.positive.definite(VarCov.test)) {
                    S.z <- try(t(chol(VarCov)), silent=TRUE)
                    if(!inherits(S.z, "try-error")) {
                         VarCov <- VarCov.test
                         S <- S.z}}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))}
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
RJ <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, bin.n, bin.p, parm.p, selectable, selected)
     {
     cur.parm <- cur.sel <- selected
     cur.parm[which(selectable == 0)] <- 1
     nonzero.post <- rep(0, LIV)
     p <- parm.p
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose a variable to include/exclude
          v.change <- sample(LIV, 1, prob=selectable)
          prop.sel <- cur.sel
          prop.parm <- cur.parm
          ### Change proposed size, but not above bin.n
          if(sum(cur.sel) < bin.n) {
               prop.sel[v.change] <- 1 - prop.sel[v.change]
               prop.parm[v.change] <- 1 - prop.parm[v.change]}
          else if(prop.sel[v.change] == 1) 
               prop.parm[v.change] <- prop.sel[v.change] <- 0
          ### Priors
          prior.cur <- sum(dbern(cur.sel, p[which(selectable == 1)], log=TRUE),
               dbinom(sum(cur.sel), bin.n, bin.p, log=TRUE))
          prior.prop <- sum(dbern(prop.sel, p[which(selectable == 1)], log=TRUE),
               dbinom(sum(prop.sel), bin.n, bin.p, log=TRUE))
          ### Hit-And-Run Proposal Parameters
          theta <- rnorm(LIV)
          theta <- theta / sqrt(sum(theta*theta))
          lambda <- runif(1)
          ### Random-Scan Componentwise Estimation (Within-Model)
          for (j in sample(which(cur.parm == 1))) {
               ### Propose new parameter values
               temp.post <- Mo0[["parm"]]
               temp.post[which(temp.post == 0)] <- nonzero.post[which(temp.post == 0)]
               temp.post[which(cur.parm == 0)] <- 0
               prop <- Mo0[["parm"]] <- temp.post
               prop[j] <- prop[j] + lambda*theta[j]
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject (Within-Model Move)
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               if(Mo0[["parm"]][j] != 0) nonzero.post[j] <- Mo0[["parm"]][j]
               Acceptance <- Acceptance + (u * (1 / sum(cur.parm)))}
          ### Random-Scan Componentwise Estimation (Between-Models)
          prop <- Mo0[["parm"]]
          prop[v.change] <- prop.sel[v.change]*(prop[v.change] +
               lambda*theta[v.change])
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
          if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
          ### Accept/Reject (Between-Models Move)
          u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]] + prior.prop -
               prior.cur)
          if(u == TRUE) {
               Mo0 <- Mo1
               cur.sel <- prop.sel
               cur.parm <- prop.parm}
          if(Mo0[["parm"]][v.change] != 0)
               nonzero.post[v.change] <- Mo0[["parm"]][v.change]
          Acceptance <- Acceptance + (u * (1 / sum(prop.parm)))
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
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
          #  VarCov <- matrix(c(5,1,1,3,1,3),2,3)
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
SAMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, parm.names, Dyn)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          if(length(staticparms) > 1) staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
Slice <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, m, w)
     {
     uni.slice <- function(x0, g, j, w=1, m=Inf, lower=-Inf, upper=Inf,
          gx0=NULL, Data)
          {
          logy <- gx0[["LP"]] - rexp(1)
          ### Find the initial interval to sample from
          u <- runif(1,0,w)
          L <- x0[j] - u
          R <- x0[j] + (w-u)
          ### Expand the interval until its ends are outside the slice,
          ### or until the limit on steps is reached
          intL <- intR <- x0
          intL[j] <- L
          intR[j] <- R
          ### Unlimited number of steps
          if(is.infinite(m)) { 
               repeat {
                    if(L <= lower) break
                    intL[j] <- L
                    gL <- g(intL, Data)
                    if(!is.finite(gL[["LP"]])) {
                         L <- L + w
                         break}
                    if(gL[["LP"]] <= logy) break
                    L <- L - w}
               repeat {
                    if(R >= upper) break
                    intR[j] <- R
                    gR <- g(intR, Data)
                    if(!is.finite(gR[["LP"]])) {
                         R <- R - w
                         break}
                    if(gR[["LP"]] <= logy) break
                    R <- R + w}
               }
          else if(m > 1) {
               ### Limited number of steps
               J <- floor(runif(1,0,m))
               K <- (m-1) - J
               while (J > 0) {
                    if(L <= lower) break
                    intL[j] <- L
                    gL <- g(intL, Data)
                    if(!is.finite(gL[["LP"]])) {
                         L <- L + w
                         break}
                    if(gL[["LP"]] <= logy) break
                    L <- L - w
                    J <- J - 1}
               while (K > 0) {
                    if(R >= upper) break
                    intR[j] <- R
                    gR <- g(intR, Data)
                    if(!is.finite(gR[["LP"]])) {
                         R <- R - w
                         break}
                    R <- R + w
                    K <- K - 1}
               }
          ### Shrink interval to lower and upper bounds
          if(L < lower) L <- lower
          if(R > upper) R <- upper
          #### Sample from the interval, shrinking it on each rejection
          prop <- x0
          repeat { 
               x1 <- runif(1,L,R)
               prop[j] <- x1
               gx1 <- g(prop, Data)
               if(is.finite(gx1[["LP"]])) {
                    x1 <- gx1[["parm"]][j]
                    if(gx1[["LP"]] >= logy) break
                    if(x1 > x0[j]) R <- x1
                    else L <- x1}
               }
          return(gx1)
          }
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Univariate Slice Sampling
               Mo0 <- Mo1 <- uni.slice(x0=Mo0[["parm"]], g=Model, j=j,
                    #w=1, m=Inf, lower=-Inf, upper=Inf, gx0=Mo0, Data)
                    w=w[j], m=m[j], lower=-Inf, upper=Inf, gx0=Mo0, Data)}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Iterations,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=apply(thinned, 2, var))
     return(out)
     }
SMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, parm.names, Dyn)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          if(length(staticparms) > 1) staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
THMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, epsilon, L, Temperature)
     {
     gr <- partial(Model, Mo0[["parm"]], Data)
     sqrt.Temp <- sqrt(Temperature)
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
          momentum1 <- momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          Mo0.1 <- Mo0
          for (l in 1:L) {
               if(2*(l-1) < L) momentum1 <- momentum1 * sqrt.Temp
               else momentum1 <- momentum1 / sqrt.Temp
               momentum1 <- momentum1 + (epsilon/2) * gr
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr <- partial(Model, prop, Data)
               momentum1 <- momentum1 + (epsilon/2) * gr
               if(2*l > L) momentum1 <- momentum1 / sqrt.Temp
               else momentum1 <- momentum1 * sqrt.Temp}
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
twalk <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, SIV, n1, at, aw)
     {
     xp0 <- SIV
     IntProd <- function(x) {return(sum(x*x))}
     DotProd <- function(x, y) {return(sum(x*y))}
     Simh1 <- function(dim, pphi, x, xp, beta)
          {
          phi <- runif(dim) < pphi
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + beta*(xp[i] - x[i]))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi)))
          }
     Simfbeta <- function(at)
          {
          if(runif(1) < (at-1)/(2*at))
               return(exp(1/(at + 1)*log(runif(1))))
          else
               return(exp(1/(1 - at)*log(runif(1))))
          }
     Simh2 <- function(dim, pphi, aw, x, xp)
          {
          u <- runif(dim)
          phi <- runif(dim) < pphi
          z <- (aw/(1+aw))*(aw*u^2 + 2*u -1)
          z <- z*phi
          return(list(rt=x + (x - xp)*z, nphi=sum(phi)))
          }
     Simh3 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))
          x + sigma*rnorm(dim)*phi
          return(list(rt=x + sigma*rnorm(dim)*phi, nphi=sum(phi),
               sigma=sigma))
          }
     G3U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd(h - xp)/(sigma^2))
          else
               return(0) 
          }
     Simh4 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))/3
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + sigma*rnorm(1))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi), sigma=sigma))
          }
     G4U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd((h - x))/(sigma^2))
          else
               return(0)
          }
     OneMove <- function(dim, Model, Data, x, U, xp, Up, at=at, aw=aw,
          pphi=pphi, F1=0.4918, F2=0.9836, F3=0.9918, Mo0.1, Mo0.2)
          {
          dir <- runif(1) ### Determine which set of points
          ker <- runif(1) ### Choose a kernel
          if(ker < F1) {
               ### Kernel h1: Traverse
               funh <- 1
               if(dir < 0.5) {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, xp, x, beta)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["LP"]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propUp <- Mo1.2[["LP"]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, x, xp, beta)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propU <- Mo1.1[["LP"]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F2) {
               ### Kernel h2: Walk
               funh <- 2
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh2(dim, pphi, aw, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh2(dim, pphi, aw, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propU <- Mo1.1[["LP"]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F3) {
               ### Kernel h3: Blow
               funh <- 3
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh3(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         W1 <- G3U(nphi, sigma,  yp, xp,  x)
                         W2 <- G3U(nphi, sigma,  xp, yp,  x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh3(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[["LP"]] * -1
                         W1 <- G3U(nphi, sigma, y, x, xp)
                         W2 <- G3U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else {
               ## Kernel h4: Hop
               funh <- 4
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh4(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         W1 <- G4U(nphi, sigma, yp, xp, x)
                         W2 <- G4U(nphi, sigma, xp, yp, x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh4(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[["LP"]] * -1
                         W1 <- G4U(nphi, sigma, y, x, xp)
                         W2 <- G4U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if(check1 & check2 & is.finite(A) & (dir < 0.5))
               Mo0.2 <- Mo1.2
          else if(check1 & check2 & is.finite(A) & (dir >= 0.5))
               Mo0.1 <- Mo1.1
          else if(!is.finite(A)) A <- 0
          #else if(!is.finite(A)) {
          #     #### Debugging
          #     cat("\ntwalk: Error in evaluating the objective:\n")
          #     cat( funh, "DU=", (U - propU), "DUp=", (Up - propUp),
          #          "DW=", (W1 - W2), "A=", A, "\n", y, "\n", yp, "\n")
          #     cat("U=", U, "propU=", propU, "Up=", Up, "propUp", propUp)}
          return(list(y=y, propU=propU, yp=yp, propUp=propUp, A=A,
               funh=funh, nphi=nphi, Mo0.1=Mo0.1, Mo0.2=Mo0.2))
          }
     Runtwalk <- function(Iterations, dim, x0, xp0, pphi, at, aw,
          F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, Model, Data, Status,
          Thinning, Acceptance, Dev, Mon, Mo0, thinned)
          {
          x <- x0 ### Primary vector of initial values
          xp <- xp0 ### Secondary vector of initial values
          Mo0.1 <- try(Model(x, Data), silent=TRUE)
          Mo0.2 <- try(Model(xp, Data), silent=TRUE)
          if(inherits(Mo0.1, "try-error") | inherits(Mo0.2, "try-error"))
               stop("Error in estimating the log-posterior.")
          if(!is.finite(Mo0.1[["LP"]]) | !is.finite(Mo0.2[["LP"]]))
               stop("The log-posterior is non-finite.")
          if(identical(x, as.vector(Mo0.1[["parm"]])) &
               identical(xp, as.vector(Mo0.2[["parm"]]))) {
               U <- Mo0.1[["LP"]] * -1
               Up <- Mo0.2[["LP"]] * -1}
          else {
               cat("\nInitial values are out of support.")
               cat("\n  Initial.Values=", x)
               cat("\n SIV=", xp)
               stop("Try re-specifying initial values.")}
          if(any(abs(x - xp) <= 0))
               stop("\nBoth vectors of initial values are not unique.")
          Acceptance <- 0
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) {cat("Iteration: ", iter,
                    ",   Proposal: Multivariate Subset\n", sep="")}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0.1[["parm"]]
                    Dev[t.iter] <- Mo0.1[["Dev"]]
                    Mon[t.iter,] <- Mo0.1[["Monitor"]]}
               ### Assign x and xp
               x <- as.vector(Mo0.1[["parm"]])
               xp <- as.vector(Mo0.2[["parm"]])
               ### Propose New Parameter Values
               move <- OneMove(dim=dim, Model, Data, x, U, xp, Up,
                    at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3,
                    Mo0.1=Mo0.1, Mo0.2=Mo0.2)
               ### Accept/Reject
               if(runif(1) < move$A) {
                    Mo0.1 <- move$Mo0.1
                    Mo0.2 <- move$Mo0.2
                    Acceptance <- Acceptance + 1 #move$nphi/dim
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- move$Mo0.1[["parm"]]
                         Dev[t.iter] <- move$Mo0.1[["Dev"]]
                         Mon[t.iter,] <- move$Mo0.1[["Monitor"]]}
                    x <- move$y
                    U <- move$propU
                    xp <- move$yp
                    Up <- move$propUp
                    }
               }
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=apply(thinned, 2, var))
          return(out)
          }
     out <- Runtwalk(Iterations=Iterations, dim=LIV, x0=Mo0[["parm"]], xp0=xp0,
          pphi=min(LIV, n1)/LIV, at=6, aw=1.5, Model=Model, Data=Data,
          Status=Status, Thinning=Thinning, Acceptance=Acceptance, Dev=Dev,
          Mon=Mon, Mo0=Mo0, thinned=thinned)
     ### Output
     return(out)
     }
USAMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, parm.names, Dyn, Fit, Begin)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.")
     ivs <- Mo0[["parm"]]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo0[["parm"]]}
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }
USMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
     thinned, tuning, parm.names, Dyn, Fit, Begin)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.")
     ivs <- Mo0[["parm"]]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[["LP"]])) Mo1 <- Mo0
               if(!is.finite(Mo1[["Dev"]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[["Monitor"]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo0[["parm"]]}
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=tuning)
     return(out)
     }

#End
