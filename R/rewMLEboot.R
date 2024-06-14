############################################################################
## rewMLE: Re-weighted MLE of BC / YJ transformation parameter lambda     ##
############################################################################

#' re-weighted Maximum Likelihood Estimator (MLE) with Bootstrapping
#' 
#' Functions rewMLE, rewMLE.boot, and rewMLE.boot.ci as custom are custom
#' wrappers derived from the "transfo" function of the cellWise package to 
#' specifically adapt the re-weighted MLE for use with clinical routine data 
#' from laboratory medicine. They contain modifications on the published
#' code from the cellWise package
#' 
#' Original Package Reference:
#' cellWise: Analyzing Data with Cellwise Outliers
#' Version: 2.5.2
#' License: GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]
#' https://cran.r-project.org/web/packages/cellWise/index.html
#' @export
rewMLE.boot = function(X, type = "YJ",
                       standardize = TRUE, 
                       quant = 0.995, nbsteps = 2,
                       winput = NULL,
                       nBoot = 200) {
  
  if (nbsteps < 1) {
    stop("nbsteps should be at least 1.")
  }
  
  if (nBoot < 1 || nBoot >= 10000){
    stop("nBoot should be an integer between 1 and 10000")
  }
  
  # Input of X:
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  X <- as.matrix(X)
  
  if (0 %in% dim(X)){
    stop("X does not contain data")
  }
  
  # Variables
  dorig       <- ncol(X) 
  colnamX     <- colnames(X)
  n           <- nrow(X)
  d           <- ncol(X)
  Xt          <- X
  
  if(is.null(colnamX)) { colnamX = seq_len(dorig) }
  
  # Setup Bootstrapping
  nBoot <- floor(nBoot)
  dimnames_list <- list(seq_len(nBoot), colnamX)
  
  # Pre-allocate Output
  lambdahats  <- matrix(1, nBoot, d, dimnames = dimnames_list)
  ttypes      <- matrix(NA, nBoot, d, dimnames = dimnames_list)
  objective   <- matrix(NA, nBoot, d, dimnames = dimnames_list)
  weights     <- matrix(0, n, d)
  nOut        <- rep(NA, d)
  
  muhat       <- matrix(0, nBoot, d, dimnames = dimnames_list)
  locx        <- matrix(0, nBoot, d, dimnames = dimnames_list)
  sigmahat    <- matrix(1, nBoot, d, dimnames = dimnames_list)
  scalex      <- matrix(1, nBoot, d, dimnames = dimnames_list)
  lambdarange <- c(-4, 4)
  
  # Setup the weights that are passed to function
  wghtsIn     <- NULL
  
  if (!is.null(winput)){
    wghtsIn <- matrix(data = winput * 1, n, d)
  }
  
  # Col-wise re-weighted MLE
  for (j in seq_len(d)) { # j=1
    x         <- X[, j]
    w         <- wghtsIn[,j]
    goodInds  <- which(!is.na(x))
    x         <- x[goodInds]
    w         <- w[goodInds]
    nOut[j]   <- length(x)
    ttype     <- getTtype(x, type, j)
    # Order
    order.x   <- order(x)
    xsort     <- x[order.x]
    wsort     <- w[order.x]
    lambdarangetemp <- lambdarange
    
    # Boostrapping
    for (i in seq_len(nBoot)) { # i=1
      ttypes[i,j]         <- ttype
      indices             <- sort(sample(order.x, size = length(order.x), replace = TRUE))
      
      if (ttype == "BC" || ttype == "bestObj") {
        lambdarangetemp   <- lambdarange
        converged         <- FALSE
        while (!converged) {
          est.out.BC      <- RewML_BC(xsort[indices], 
                                      lambdarange = lambdarangetemp, 
                                      standardize = standardize, 
                                      init = "BCr", 
                                      winput = wsort[indices],
                                      quant = quant, nbsteps = nbsteps)
          
          converged       <- min(abs(est.out.BC$lambdahat.rew - 
                                       lambdarangetemp)) > diff(lambdarangetemp) * 
            0.05
          
          lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
            (1 + (abs(est.out.BC$lambdahat.rew - lambdarangetemp) == 
                    min(abs(est.out.BC$lambdahat.rew - lambdarangetemp))) + 
               0)
        }
        lambdahats[i,j]               <- est.out.BC$lambdahat.rew
        objective[i,j]                <- est.out.BC$critval.rew
        Xt[goodInds[indices], j]      <- est.out.BC$yt.rew
        weights[goodInds[indices], j] <- est.out.BC$weights + weights[goodInds[indices], j]
        muhat[i,j]                    <- est.out.BC$muhat
        sigmahat[i,j]                 <- est.out.BC$sigmahat
        locx[i,j]                     <- est.out.BC$locx
        scalex[i,j]                   <- est.out.BC$scalex
      }
      if (ttype == "YJ" || ttype == "bestObj") {
        lambdarangetemp     <- lambdarange
        converged           <- FALSE
        while (!converged) {
          est.out.YJ        <- RewML_YJ(xsort[indices],
                                        lambdarange = lambdarangetemp, 
                                        standardize = standardize, 
                                        init = "YJr", 
                                        winput = wsort[indices],
                                        quant = quant, nbsteps = nbsteps)
          
          converged         <- min(abs(est.out.YJ$lambdahat.rew - 
                                         lambdarangetemp)) > diff(lambdarangetemp) * 
            0.05
          lambdarangetemp   <- 1 + (lambdarangetemp - 1) * 
            (1 + (abs(est.out.YJ$lambdahat.rew - lambdarangetemp) == 
                    min(abs(est.out.YJ$lambdahat.rew - lambdarangetemp))) + 
               0)
        }
        lambdahats[i,j]                 <- est.out.YJ$lambdahat.rew
        objective[i,j]                  <- est.out.YJ$critval.rew
        Xt[goodInds[indices], j]        <- est.out.YJ$yt.rew
        weights[goodInds[indices], j]   <- est.out.YJ$weights + weights[goodInds[indices], j]
        muhat[i,j]                      <- est.out.YJ$muhat
        sigmahat[i,j]                   <- est.out.YJ$sigmahat
        locx[i,j]                       <- est.out.YJ$locx
        scalex[i,j]                     <- est.out.YJ$scalex
      }
      if (ttype == "bestObj") {
        if (est.out.BC$critval.rew < est.out.YJ$critval.rew) {
          lambdahats[i,j]               <- est.out.BC$lambdahat.rew
          objective[i,j]                <- est.out.BC$critval.rew
          Xt[goodInds[order.x], j]      <- est.out.BC$yt.rew
          weights[goodInds[indices], j] <- est.out.BC$weights + weights[goodInds[indices], j]
          muhat[i,j]                    <- est.out.BC$muhat
          sigmahat[i,j]                 <- est.out.BC$sigmahat
          locx[i,j]                     <- est.out.BC$locx
          scalex[i,j]                   <- est.out.BC$scalex
          ttypes[i,j]                   <- "BC"
        }
        else {
          ttypes[i,j]                   <- "YJ"
        }
      }
    }
    
  }
  Y <- if (standardize) { scale(Xt, center = colMeans(muhat), scale = colMeans(sigmahat))
  } else { Xt }
  
  # Standardize weights
  weights <-    weights/max(weights, na.rm = T)
  wthresh <-    rep(0, d)
  wbinary <-    as.matrix(weights)
  
  for (j in seq_len(d)){
    wcol <-     weights[,j]
    wthresh[j] <-  mean(wcol[wcol > 0], na.rm=T) - 1.5*sd(wcol[wcol > 0], na.rm=T)
    wbinary[,j]   <- ifelse(wcol > wthresh[j], 1, 0)
    
  }
  
  # Results
  result                <- c(list(nBoot = nBoot,
                                  lambdahats = lambdahats, 
                                  objective = objective, 
                                  Y = Y, X = X, Xt = Xt,
                                  weights = weights,
                                  wbinary = wbinary,
                                  wthresh = wthresh,
                                  
                                  nOut = nOut,
                                  ttypes = ttypes, 
                                  muhat = muhat, 
                                  sigmahat = sigmahat, 
                                  locx = locx, 
                                  scalex = scalex, 
                                  standardize = standardize),
                             list(dorig = dorig, colnamX = colnamX,
                                  quant = quant, nbsteps = nbsteps,
                                  winput = winput) )
  
  # Assigning the rewMLE.boot class
  class(result)   <- "rewMLE.boot"
  return(result)
  
}


#' Estimate the adjusted bootstrap percentile (BCa)
#' 
#' Calculating the 90% adjusted bootstrap percentile (BCa). If the estimator
#' is sensitive to the data during bootstrapping, this fixes the skewness 
#' and the apparent outliers of a distribution of bootstrapped estimates
#' @export
estBCa <- function(boot.est, orig_stat, CI = 0.9) {
  
  # Calculate Bias Correction (z_0)
  p_hat <- mean(boot.est < orig_stat)
  z_0 <- qnorm(p_hat)
  
  # Jackknife estimator (function estJK)
  jk_result <- estJK(boot.est)
  a         <- jk_result$acceleration
  
  # Calculate Adjusted Percentiles for BCa CI
  alpha     <- 1-CI # For 90% CI
  z_alpha1  <- qnorm(alpha / 2)
  z_alpha2  <- qnorm(1 - alpha / 2)
  p1        <- pnorm(z_0 + (z_0 + z_alpha1) / (1 - a * (z_0 + z_alpha1)))
  p2        <- pnorm(z_0 + (z_0 + z_alpha2) / (1 - a * (z_0 + z_alpha2)))
  
  lower_bca <- quantile(boot.est, p1)
  upper_bca <- quantile(boot.est, p2)
  
  return(c(lower_bca, upper_bca))
}


#' Jackknife Estimator
#' 
estJK <- function(data) {
  
  # Compute the original statistic for the full dataset
  orig_stat <- mean(data)
  
  # Compute jackknife estimates
  n <- length(data) 
  jkhats <- numeric(n) 
  
  for (i in 1:n) {
    
    # Excluding the i-th observation and compute
    # the statistic for this subset
    subset_data <- data[-i]
    jkhats[i] <- mean(subset_data)
    
  }
  
  # Compute the mean of the jackknife estimates
  jkmean <- mean(jkhats)
  
  # Compute the acceleration (a)
  num <- sum((jkmean - jkhats)^3)
  dom <- 6 * (sum((jkmean - jkhats)^2)^(3/2))
  a <- num / dom
  
  return(list(estJK = jkhats, acceleration = a))
}

#' Confidence Intervals for the rewMLE Bootstrapping
#' 
#' @export
rewMLE.boot.ci = function(rewMLE.boot.out, CI = 0.9){
  
  lambdahats  <- rewMLE.boot.out$lambdahats
  names       <- colnames(lambdahats)
  nCols       <- ncol(lambdahats)
  alpha       <- 1 - CI
  
  muhats      <- rewMLE.boot.out$muhat
  sigmahats   <- rewMLE.boot.out$sigmahat
  locxs       <- rewMLE.boot.out$locx
  scalexs     <- rewMLE.boot.out$scalex
  
  
  pass        <- list("objective" = rewMLE.boot.out$objective, 
                      "Y" = rewMLE.boot.out$Y, "X" = rewMLE.boot.out$X,
                      "ttypes" = rewMLE.boot.out$ttypes, 
                      "standardize" = rewMLE.boot.out$standardize)
  
  # Weights
  weights     <- rewMLE.boot.out$weights
  
  
  # Calculate the CIs using BCa
  lambdaEst   <- matrix(0, nrow=2, ncol=nCols,
                        dimnames = list(c("mean", "median"), names))
  lambdaBCa   <- matrix(0, nrow=2, ncol=nCols,
                        dimnames = list(c(alpha/2, 1-(alpha/2)), names))
  lambdaCIs   <- matrix(0, nrow=2, ncol=nCols,
                        dimnames = list(c(alpha/2, 1-(alpha/2)), names))
  
  muEst       <- rep(0, nCols)
  sigmaEst    <- rep(0, nCols)
  locxEst     <- rep(0, nCols)
  scalex      <- rep(0, nCols)
  wThresh     <- rep(0, nCols)
  
  for (j in seq_len(nCols)){ #j
    
    # Given the weights are cont. between 0,1
    # Round it back to binary weights based on
    wBoot         <- weights[,j]
    wThresh[j]    <- mean(wBoot[wBoot > 0], na.rm=T) - 1.5*sd(wBoot[wBoot > 0], na.rm=T)
    
    
    weights[,j]   <- ifelse(wBoot > wThresh[j], 1, 0)
    
    # Bootstrapped lambdas
    lambdasBoot   <- lambdahats[,j]
    # Calculated Mean & Median
    mean          <- mean(lambdasBoot, na.rm=T)
    median        <- median(lambdasBoot, na.rm=T)
    lambdaEst[,j] <- c(mean, median)
    
    # Compute BCa CI for this column
    BCa           <- estBCa(lambdasBoot, mean, CI = CI)
    lambdaCIs[,j] <- quantile(lambdasBoot, p = c(alpha/2, 1-(alpha/2)), type=8)
    lambdaBCa[,j] <- BCa
    
    # Decide which of the values to keep
    keep          <- which( lambdasBoot > BCa[1] & lambdasBoot <= BCa[2] )
    
    muEst[j]      <- median(muhats[,j][keep], na.rm=TRUE)
    sigmaEst[j]   <- median(sigmahats[,j][keep], na.rm=TRUE)
    locxEst[j]    <- median(locxs[,j][keep], na.rm=TRUE)
    scalex[j]     <- median(scalexs[,j][keep], na.rm=TRUE)
    
  }
  boot.ci <- c(list(lambdahats = lambdaEst,
                    lambdaCIs = lambdaCIs,
                    lambdaBCa = lambdaBCa,
                    weights = weights,
                    wThresh = wThresh,
                    muhat = muEst, 
                    sigmahat = sigmaEst, 
                    locx = locxEst, 
                    scalex = scalex),
               pass)
  
  # Assigning the rewML class
  class(boot.ci) <- "rewMLE.boot.ci"
  return(boot.ci)
}



# --------------------------------
#' print method for rewMLE.boot.ci
#' @export
print.rewMLE.boot  <- function(x, CI=0.9, ...) {
  # Print methods do not need to use '...',
  # but it's available if needed.

  # Print output
  cat("print class for `rewMLE.boot`\n\n")
  
  cat(str(x))

  invisible(x)
}

# -------------------------------
#' plot method for rewMLE.boot.ci
#' @export
plot.rewMLE.boot <- function(x, CI=0.9, trim = 0.95, ... , main = NULL){
  
  x.ci <- rewMLE.boot.ci(rewMLE.boot.out = x, CI = CI)
  nBoot <- x$nBoot
  
  weights <- x$weights
  nCols <- ncol(x$weights)
  nRows <- nrow(x$weights)
  wThresh <-  x.ci$wThresh
  
  for (i in 1:nCols){
    
    plot(x= 1:nRows, y=weights[,i],main = "Final weights")
    abline(h = wThresh[i], lty = 2)
    
  }

}
