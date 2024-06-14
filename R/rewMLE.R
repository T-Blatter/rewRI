############################################################################
## rewMLE: Re-weighted MLE of BC / YJ transformation parameter lambda     ##
############################################################################

#' re-weighted Maximum Likelihood Estimator (MLE)
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
rewMLE = function(X, type = "YJ",
                  standardize = TRUE, 
                  quant = 0.995, nbsteps = 2,
                  winput = NULL) {
  
  if (nbsteps < 1) {
    stop("nbsteps should be at least 1.")
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
  
  if (is.null(colnamX)) { colnamX = seq_len(dorig) }
  
  # Pre-allocate Output
  lambdahats  <- rep(1, d)
  ttypes      <- rep(NA, d)
  objective   <- rep(NA, d)
  weights     <- matrix(0, n, d)
  nOut        <- rep(NA, d)
  
  muhat       <- rep(0, d)
  locx        <- rep(0, d)
  sigmahat    <- rep(1, d)
  scalex      <- rep(1, d)
  lambdarange <- c(-4, 4)
  
  # Setup the weights that are passed to function
  wghtsIn     <- NULL
  
  if (!is.null(winput)){
    wghtsIn   <- matrix(data = winput * 1, n, d)
  }
  
  # ColWise re-weighted MLE
  for (j in seq_len(d)) { # j=1
    x         <- X[, j]
    w         <- wghtsIn[,j]
    goodInds  <- which(!is.na(x))
    x         <- x[goodInds]
    w         <- w[goodInds]
    nOut[j]   <- length(x)
    ttype     <- getTtype(x, type, j)
    ttypes[j] <- ttype
    # Order
    order.x   <- order(x)
    xsort     <- x[order.x]
    wsort     <- w[order.x]
    
    if (ttype == "BC" || ttype == "bestObj") {
      lambdarangetemp <- lambdarange
      converged       <- FALSE
      while (!converged) {
        est.out.BC    <- RewML_BC(xsort, 
                                  lambdarange = lambdarangetemp, 
                                  standardize = standardize, 
                                  init = "BCr", 
                                  winput = wsort,
                                  quant = quant, nbsteps = nbsteps)
        
        converged     <- min(abs(est.out.BC$lambdahat.rew - 
                                   lambdarangetemp)) > diff(lambdarangetemp) * 
          0.05
        
        lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
          (1 + (abs(est.out.BC$lambdahat.rew - lambdarangetemp) == 
                  min(abs(est.out.BC$lambdahat.rew - lambdarangetemp))) + 
             0)
      }
      lambdahats[j]                 <- est.out.BC$lambdahat.rew
      objective[j]                  <- est.out.BC$critval.rew
      Xt[goodInds[order.x], j]      <- est.out.BC$yt.rew
      weights[goodInds[order.x], j] <- est.out.BC$weights
      muhat[j]                      <- est.out.BC$muhat
      sigmahat[j]                   <- est.out.BC$sigmahat
      locx[j]                       <- est.out.BC$locx
      scalex[j]                     <- est.out.BC$scalex
    }
    if (ttype == "YJ" || ttype == "bestObj") {
      lambdarangetemp   <- lambdarange
      converged         <- FALSE
      while (!converged) {
        est.out.YJ      <- RewML_YJ(xsort,
                                    lambdarange = lambdarangetemp, 
                                    standardize = standardize, 
                                    init = "YJr", 
                                    winput = wsort,
                                    quant = quant, nbsteps = nbsteps)
        
        converged       <- min(abs(est.out.YJ$lambdahat.rew - 
                                     lambdarangetemp)) > diff(lambdarangetemp) * 
          0.05
        lambdarangetemp <- 1 + (lambdarangetemp - 1) * 
          (1 + (abs(est.out.YJ$lambdahat.rew - lambdarangetemp) == 
                  min(abs(est.out.YJ$lambdahat.rew - lambdarangetemp))) + 
             0)
      }
      lambdahats[j]                   <- est.out.YJ$lambdahat.rew
      objective[j]                    <- est.out.YJ$critval.rew
      Xt[goodInds[order.x], j]        <- est.out.YJ$yt.rew
      weights[goodInds[order.x], j]   <- est.out.YJ$weights
      muhat[j]                        <- est.out.YJ$muhat
      sigmahat[j]                     <- est.out.YJ$sigmahat
      locx[j]                         <- est.out.YJ$locx
      scalex[j]                       <- est.out.YJ$scalex
    }
    if (ttype == "bestObj") {
      if (est.out.BC$critval.rew < est.out.YJ$critval.rew) {
        lambdahats[j]                 <- est.out.BC$lambdahat.rew
        objective[j]                  <- est.out.BC$critval.rew
        Xt[goodInds[order.x], j]      <- est.out.BC$yt.rew
        weights[goodInds[order.x], j] <- est.out.BC$weights
        muhat[j]                      <- est.out.BC$muhat
        sigmahat[j]                   <- est.out.BC$sigmahat
        locx[j]                       <- est.out.BC$locx
        scalex[j]                     <- est.out.BC$scalex
        ttypes[j]                     <- "BC"
      }
      else {
        ttypes[j]                     <- "YJ"
      }
    }
  }
  Y <- if (standardize) { scale(Xt, center = muhat, scale = sigmahat)}
  else { Xt }
  
  # Results
  result               <- c(list(lambdahats = lambdahats, 
                                 objective = objective, 
                                 X = X, Y = Y, Xt = Xt,
                                 weights = weights,
                                 nOut = nOut,
                                 ttypes = ttypes, 
                                 muhat = muhat, 
                                 sigmahat = sigmahat, 
                                 locx = locx, 
                                 scalex = scalex, 
                                 standardize = standardize), 
                            list(dorig = dorig,colnamX = colnamX,
                                 quant = quant, nbsteps = nbsteps,
                                 winput = winput))
  
  # Assigning the rewMLE class
  class(result)   <- "rewMLE"
  return(result)
  
}

#' Transform data by re-weighted MLE of BC/YJ
#' 
#' This function takes new data (Xnew) and transform them
#' according to the parameters in the previously established
#' re-weighted MLE (rewMLE.out).
#' The rewMLE.out needs to an objects of class 
#' `rewMLE`, `rewMLE.boot`, or `rewMLE.boot.ci`.
#' Source: This is a modification on the `transfo_newdata` function.
#' cellWise: Analyzing Data with Cellwise Outliers
#' Version: 2.5.2
#' https://cran.r-project.org/web/packages/cellWise/index.html
transfo_rewMLE = function(Xnew, rewMLE.out) {
  if (is.vector(Xnew)) {
    Xnew <- matrix(Xnew, ncol = 1)
  }
  Xnew <- as.matrix(Xnew)
  d = ncol(Xnew)
  dorig = rewMLE.out$dorig
  # if(d != dorig) stop(paste0(
  #   "Xnew should have ",dorig," columns like the original data."))
  colnamnew = colnames(Xnew)
  if(is.null(colnamnew)) { colnamnew = seq_len(d) }
  # if(all.equal(colnamnew, rewMLE.out$colnamX) != TRUE) stop(
  #   "Xnew should have the same column names as the original data.")
  colInA <- seq_len(d)
  if (!is.null(rewMLE.out$out$colInAnalysis)) {
    colInA <- rewMLE.out$out$colInAnalysis
  }
  standardize <- rewMLE.out$standardize
  locx <- rewMLE.out$locx
  scalex <- rewMLE.out$scalex
  muhat <- rewMLE.out$muhat
  sigmahat <- rewMLE.out$sigmahat
  ttypes <- rewMLE.out$ttypes
  lambdas <- rewMLE.out$lambdahats
  Ynew <- Xnew
  for (j in colInA) {
    xnewt <- Xnew[, colInA[j]] # actual column
    xnewt <- xnewt[!is.na(xnewt)]
    if (ttypes[j] == "YJ") {
      ynewt <- xnewt
      if (standardize) {
        ynewt <- scale(ynewt, locx[j], scalex[j])
      }
      ynewt <- YJ(ynewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        ynewt <- scale(ynewt, muhat[j], sigmahat[j])
      }
      Ynew[1:length(ynewt), j] <- ynewt
    }
    else if (ttypes[j] == "BC") {
      ynewt <- xnewt[xnewt>0]
      if (standardize) {
        ynewt <- ynewt/scalex[j]
      }
      if(min(ynewt) <= 0) stop(paste0(
        "Column ",j," has some value(s) <= 0, but in the original data\n",
        "it was transformed by Box-Cox. You have to either make all\n",
        "values in this column strictly positive, or apply transfo()\n",
        "again to the original data but with type = \"YJ\"."))
      ynewt <- BC(ynewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        ynewt <- scale(ynewt, muhat[j], sigmahat[j])
      }
      Ynew[1:length(ynewt), j] <- ynewt
    }
    else {
      stop(paste0("Invalid transformation type ",ttypes[j],
                  " for column ",j))
    }
  }
  return(Ynew)
}

#' Back transform data by re-weighted MLE of BC/YJ
#' 
#' This function back-transforms data (Ynew) and transform them
#' back to the original space. This is done in accordance to the parameters 
#' in the previously established re-weighted MLE (rewMLE.out).
#' 
#' The rewMLE.out needs to an objects of class 
#' `rewMLE`, `rewMLE.boot`, or `rewMLE.boot.ci`.
#' Source: This is a modification on the `transfo_transformback` function.
#' cellWise: Analyzing Data with Cellwise Outliers
#' Version: 2.5.2
#' https://cran.r-project.org/web/packages/cellWise/index.html
transfoback_rewMLE = function (Ynew, rewMLE.out) {
  if (is.vector(Ynew)) {
    Ynew <- matrix(Ynew, ncol = 1)
  }
  Ynew <- as.matrix(Ynew)
  d = ncol(Ynew)
  dorig = rewMLE.out$dorig
  # if(d != rewMLE.out$dorig) stop(paste0(
  #   "Ynew should have ",dorig," columns like the original ",
  #   "transformed data."))
  colnamnew = colnames(Ynew)
  if(is.null(colnamnew)) { colnamnew = seq_len(d) }
  # if(all.equal(colnamnew, rewMLE.out$colnamX) != TRUE) stop(paste0(
  #   "Ynew should have the same column names as the original ",
  #   "transformed data."))
  standardize <- rewMLE.out$standardize
  locx <- rewMLE.out$locx
  scalex <- rewMLE.out$scalex
  ttypes <- rewMLE.out$ttypes
  lambdas <- rewMLE.out$lambdahats
  muhat <- rewMLE.out$muhat
  sigmahat <- rewMLE.out$sigmahat
  colInAnalysis <- seq_len(ncol(Ynew))
  # if (!is.null(rewMLE.out$out$colInAnalysis)) {
  #   colInAnalysis <- rewMLE.out$out$colInAnalysis
  # }
  Xnew <- Ynew # initialization
  for (j in colInAnalysis) {
    # Shouldn't we step over j in seq_len(length(colInAnalysis)) ?
    ynewt <- Ynew[, j]
    ynewt <- ynewt[!is.na(ynewt)]
    if (ttypes[j] == "YJ") {
      xnewt <- ynewt
      stan = " "
      if (standardize) {
        xnewt <- xnewt * sigmahat[j] + muhat[j]
        stan = " destandardized "
      }
      if(lambdas[j] > 2){
        lowerb = -1/abs(lambdas[j] - 2)
        newlowerb = 0.95*lowerb
        if(min(xnewt) < lowerb) warning(paste0(
          "The lowest",stan,"value in column ",j," was ",min(xnewt),
          " .\n","This is below the expected lower bound ",lowerb,
          " .\n","Such low",stan,"values were put equal to ",
          newlowerb,"\nso they could be transformed back."))
        xnewt = pmax(xnewt, newlowerb)
      }
      if(lambdas[j] < 0){
        upperb = 1/abs(lambdas[j])
        newupperb = 0.95*upperb
        if(max(xnewt) > upperb) warning(paste0(
          "The highest",stan,"value in column ",j," was ",max(xnewt),
          " .\n","This is above the expected upper bound ",upperb,
          " .\n","Such high",stan,"values were put equal to ",
          newupperb,"\nso they could be transformed back."))
        xnewt = pmin(xnewt, newupperb)
      }
      xnewt <- iYJ(xnewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        xnewt <- xnewt * scalex[j] + locx[j]
      }
      Xnew[1:length(xnewt), j] <- xnewt
    }
    else if (ttypes[j] == "BC") {
      xnewt <- ynewt
      stan = " "
      if (standardize) {
        xnewt <- xnewt * sigmahat[j] + muhat[j]
        stan = " destandardized "
      }
      if(lambdas[j] > 0){
        lowerb = -1/abs(lambdas[j])
        newlowerb = 0.95*lowerb
        if(min(xnewt) < lowerb) warning(paste0(
          "The lowest",stan,"value in column ",j," was ",min(xnewt),
          " .\n","This is below the expected lower bound ",lowerb,
          " .\n","Such low",stan,"values were put equal to ",
          newlowerb,"\nso they could be transformed back."))
        xnewt = pmax(xnewt, newlowerb)  
      }
      if(lambdas[j] < 0){
        upperb = 1/abs(lambdas[j])
        newupperb = 0.95*upperb
        if(max(xnewt) > upperb) warning(paste0(
          "The highest",stan,"value in column ",j," was ",max(xnewt),
          " .\n","This is above the expected upper bound ",upperb,
          " .\n","Such high",stan,"values were put equal to ",
          newupperb,"\nso they could be transformed back."))
        xnewt = pmin(xnewt, newupperb)
      }
      xnewt <- iBC(xnewt, lambdas[j], stdToo = FALSE)$yt
      if (standardize) {
        xnewt <- xnewt * scalex[j]
      }
      Xnew[1:length(xnewt), j] <- xnewt
    }
    else {
      stop(paste0("Invalid transformation type ",ttypes[j],
                  " for column ",j))
    }
  }
  return(Xnew)
}

#' The reweighted ML estimator for the Box-Cox transformation parameter.
#'
#' args: 
#'  x:              vector of _sorted_ observations
#'  lambdarange:    grid of lambda values. If NULL, a grid between
#'                  -2 and 4 is chosen
#'  init:           initial estimator. should be "BCr" or "BC"
#'  quant:         quantile for determining the weights in the 
#'                  reweighting step
#'  nbsteps:        number of reweighting steps
RewML_BC <- function(x,
                     lambdarange = NULL,
                     standardize = TRUE,
                     init = "BCr",
                     winput = NULL,
                     quant = 0.99,
                     nbsteps = 2) {
  
  if (!init %in% c("BC", "BCr")) {
    stop("init should be either 'BC' or 'BCr'")
  }
  
  x <- na.omit(x)
  
  if(!is.null(winput)){
    w <- winput
    w[is.na(winput)] <- 0
  } else { w <- NULL}
  
  
  if (is.unsorted(x)) {
    x.order <- order(x, decreasing = FALSE)
    x       <- x[x.order]
    w       <- w[x.order]
  } else {
    x.order <- NULL
  }
  
  scalex <- 1
  if (standardize) {
    scalex <- median(x)
    x <- x / scalex
  }
  
  # Range of lambda over which to optimize:
  if (is.null(lambdarange)) { lambdarange <- c(-4, 6) }
  if (init == "BCr") {
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = BCr)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    BCr.out.raw   <- BCr(x, lambdahat.raw)
    yt.raw        <- BCr.out.raw$yt
    zt.raw        <- BCr.out.raw$zt
  } else { 
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = BC)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    BC.out.raw    <- BC(x, lambdahat.raw)
    yt.raw        <- BC.out.raw$yt
    zt.raw        <- BC.out.raw$zt
  }
  
  rew.out <- reweightBCr(x = x, zt.raw = zt.raw,
                         lambdahat.raw = lambdahat.raw, 
                         lambdarange = lambdarange, quant = quant,
                         winput = w,
                         nbsteps = nbsteps,
                         standardize = standardize)
  
  BC.out.rew    <- rew.out$BC.out.rew
  yt.rew        <- BC.out.rew$yt
  zt.rew        <- BC.out.rew$zt
  critval.rew   <- rew.out$critval.rew
  lambdahat.rew <- rew.out$lambdahat.rew
  weights       <- rew.out$wgts
  muhat         <- mean(yt.rew[weights == 1])
  sigmahat      <- sd(yt.rew[weights == 1])
  
  if (!is.null(x.order)) {
    yt.raw[x.order] <- yt.raw
    zt.raw[x.order] <- zt.raw
    yt.rew[x.order] <- yt.rew
    zt.rew[x.order] <- zt.rew
    weights[x.order] <- weights
  }
  return(list(lambdahat.raw = lambdahat.raw,
              yt.raw = yt.raw,
              zt.raw = zt.raw,
              weights = weights,
              lambdahat.rew = lambdahat.rew,
              yt.rew = yt.rew,
              zt.rew = zt.rew,
              critval.rew = critval.rew, 
              muhat = muhat,
              sigmahat = sigmahat,
              locx = 0,
              scalex = scalex))
}

#' The reweighted ML estimator for the Yeo-Johnson transformation parameter.
#'
#' args: 
#'   x             : vector of _sorted_ observations
#'   lambdarange       : grid of lambda values. If NULL, a grid between 
#'                   -2 and 4 is chosen
#'   init          : initial estimator. should be "YJr" or "YJ"
#'   quant         : quantile for determining the weights in the
#'                   reweighting step
#'   nbsteps       : number of reweighting steps
RewML_YJ <- function(x, 
                     lambdarange = NULL, 
                     quant = 0.99, 
                     init = "YJr",
                     winput = NULL,
                     standardize = TRUE, 
                     nbsteps = 2) {
  
  if (!init %in% c("YJ", "YJr")) {
    stop("init should be either 'YJ' or 'YJr'")
  }
  
  x <- na.omit(x)
  
  if(!is.null(winput)){
    w <- winput
    w[is.na(winput)] <- 0
  } else {
    w <- NULL
  }
  
  
  if(is.unsorted(x)) {
    x.order <- order(x, decreasing = FALSE)
    x       <- x[x.order]
    w       <- w[x.order]
  } else {
    x.order <- NULL
  }
  
  locx <- 0
  scalex <- 1
  if (standardize) {
    locx <- median(x)
    scalex <- mad(x)
    x <- (x - locx) / scalex
  }
  
  # Range of lambda over which to optimize:
  if (is.null(lambdarange)) { lambdarange <- c(-4,6) }
  if (init == "YJr") {
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = YJr)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    YJr.out.raw   <- YJr(x, lambdahat.raw)
    yt.raw        <- YJr.out.raw$yt
    zt.raw        <- YJr.out.raw$zt
  } else { 
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = YJ)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    YJr.out.raw   <- YJ(x, lambdahat.raw)
    yt.raw        <- YJr.out.raw$yt
    zt.raw        <- YJr.out.raw$zt
  }
  # Reweighting:
  rew.out <- reweightYJr(x, zt.raw, lambdahat.raw, lambdarange,
                         winput = w, 
                         quant = quant, nbsteps = nbsteps)
  weights       <- rew.out$wgts
  YJ.out.rew    <- rew.out$YJ.out.rew
  yt.rew        <- YJ.out.rew$yt
  zt.rew        <- YJ.out.rew$zt
  critval.rew   <- rew.out$critval.rew
  lambdahat.rew <- rew.out$lambdahat.rew
  muhat         <- mean(yt.rew[weights == 1])
  sigmahat      <- sd(yt.rew[weights == 1])
  
  if (!is.null(x.order)) {
    yt.raw[x.order] <- yt.raw
    zt.raw[x.order] <- zt.raw
    yt.rew[x.order] <- yt.rew
    zt.rew[x.order] <- zt.rew
    weights[x.order] <- weights
  }
  
  return(list(lambdarange = lambdarange,
              lambdahat.raw = lambdahat.raw,
              yt.raw = yt.raw,
              zt.raw = zt.raw,
              weights = weights,
              critval.rew = critval.rew,
              lambdahat.rew = lambdahat.rew,
              yt.rew = yt.rew,
              zt.rew = zt.rew, 
              muhat = muhat,
              sigmahat = sigmahat,
              locx = locx,
              scalex = scalex))
}

#' Function for reweighted maximum likelihood, based on an initial estimate.
#' args: 
#'   x             : vector of sorted original observations
#'   zt.raw        : vector of sorted poststandardized transformed data 
#'                   (with initial lambdahat)
#'   lambdahat.raw : initial estimate for transformation parameter lambda
#'   lambdarange   : range of lambda values.
#'   quant         : quantile for determining the weights in the 
#'                   reweighting step
#'   nbsteps       : number of reweighting steps
reweightYJr <- function(x, 
                        zt.raw, 
                        lambdahat.raw,
                        lambdarange,
                        winput = NULL, 
                        quant = 0.99, 
                        nbsteps = 2) {
  
  if (is.null(lambdarange)) { lambdarange = c(-4,6) }
  
  # Initial re-weighting:
  # TODO: Check, whether this is done correctly
  if(is.null(winput) || sum(winput) == 0){
    weights  <- abs(zt.raw) <= sqrt(qchisq(quant,1))
  } else {
    weights  <- (winput*1) ##| abs(zt.raw) <= sqrt(qchisq(quant,1))
  }
  
  rewinds <- which(weights == 1)
  x.rew   <- x[rewinds]
  x.all   <- x
  
  # YJ transform with weighted ML:
  for (k in seq_len(nbsteps)) {
    lambdahat.rew <- estML(x = x.rew, lambdarange = lambdarange, 
                           type = "YJ", standardize = FALSE)$lambda
    YJ.out.rew    <- YJ(x.all, lambdahat.rew)
    wgts    <- abs(YJ.out.rew$zt) <= sqrt(qchisq(quant,1))
    rewinds <- which(wgts == 1)
    x.rew   <- x[rewinds]
  }
  critval.rew <- getCritval(x.all, lambdahat.rew, tfunc = YJ,
                                       quant=quant)
  return(list(wgts = wgts,
              critval.rew  = critval.rew,
              lambdahat.rew = lambdahat.rew,
              YJ.out.rew = YJ.out.rew))
}

#' function for reweighted maximum likelihood, based on an initial estimate.
#' args: 
#'   x              : vector of sorted original observations
#'   zt.raw         : vector of sorted poststandardized transformed data 
#'                    (from the initial lambdahat)
#'   lambdahat.raw  : initial estimate for transformation parameter lambda
#'   lambdarange    : range of lambda values.
#'   quant          : quantile for determining the weights in the 
#'                    reweighting step
#'   nbsteps        : number of reweighting steps
reweightBCr <- function(x, 
                        zt.raw, 
                        lambdahat.raw, 
                        lambdarange, 
                        winput = NULL,
                        quant = 0.99,
                        nbsteps = 2, standardize = TRUE) {
  
  if (is.null(lambdarange)) { lambdarange = c(-4,6) }
  
  if(is.null(winput) || sum(winput) == 0){
    weights  <- abs(zt.raw) <= sqrt(qchisq(quant,1))
  } else {
    weights  <- (winput*1) ##| abs(zt.raw) <= sqrt(qchisq(quant,1))
  }
  
  x.all   <- x
  x.rew   <- x[which(weights == 1)]
  
  # BC transform with weighted ML
  for (k in seq_len(nbsteps)) {
    lambdahat.rew <- estML(x = x.rew, lambdarange = lambdarange, 
                           type = "BC", standardize = FALSE)$lambda
    BC.out.rew <- BC(x.all, lambdahat.rew)
    wgts      <- abs(BC.out.rew$zt) <= sqrt(qchisq(quant,1))
    rewinds    <- which(wgts == 1)
    x.rew      <- x[rewinds]
  }
  critval.rew <- getCritval(x.all, lambdahat.rew, tfunc = BC,quant = quant)
  return(list(BC.out.rew = BC.out.rew,
              lambdahat.rew = lambdahat.rew, 
              critval.rew = critval.rew,
              wgts = wgts))
}

#' The Box-Cox transformation.
BC <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  
  if (lambda == 0) {
    yt <- log(y)
  } else {
    yt <- (y^lambda - 1) / lambda
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(yt, ncol = 1), 
                                        type = "hubhub")
      zt <- (yt - locScale$loc) / locScale$scale
    } else { 
      zt <- yt
    }
  } else {
    zt <- NULL
  }
  return(list(yt = yt, zt =zt))
}

#' Inverse of the Box-Cox transformation.
iBC <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  if (lambda == 0) {
    yt <- exp(y)
  } else {
    yt <- (1 + y * lambda)^(1/lambda)
  }
  
  if (stdToo) {
    zt = (yt - median(yt)) / mad(yt)
  } else {
    zt = NULL
  }
  return(list(yt = yt, zt = zt))
}

#' The Yeo-Johnson transformation.
YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  indlow  <- which(y < 0)
  indhigh <- which(y >= 0)
  if (lambda != 0) {
    y[indhigh] = ((1 + y[indhigh])^(lambda) - 1) / lambda
  } else { 
    y[indhigh] = log(1 + y[indhigh])
  }
  if (lambda != 2) {
    y[indlow] = -((1 - y[indlow])^(2 - lambda) - 1) / (2 - lambda)
  } else {
    y[indlow] = -log(1 - y[indlow])
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(y, ncol = 1), 
                                        type = "hubhub")
      zt <- (y - locScale$loc) / locScale$scale
    } else {
      zt <- y
    }
  } else {
    zt <- NULL
  }
  return(list(yt = y, zt = zt))
}

#' Inverse of the Yeo-Johnson transformation.
iYJ <- function(y, lambda, stdToo = TRUE) {
  indlow  <- which(y < 0)
  indhigh <- which(y >= 0)
  if (lambda != 0) {
    y[indhigh] = (1 + lambda * y[indhigh])^(1 / lambda) - 1
  } else { 
    y[indhigh] = exp(y[indhigh]) - 1
  }
  if (lambda != 2) {
    y[indlow] = -((1 + (lambda - 2) * y[indlow])^(1/(2-lambda))) + 1
  } else {
    y[indlow] = 1 - exp(-y[indlow])
  }
  if (stdToo) {
    zt = (y - median(y)) / mad(y)
  } else {
    zt = NULL
  }
  return(list(yt = y, zt = zt))
}

#' partial derivative of YJ with respect to the argument x.
YJprimex <- function(x, lambda) {
  (1 + abs(x))^(sign(x) * (lambda - 1))
}


BCr = function(y, lambda, chg = NULL, stdToo = TRUE) {
  # Rectified Box Cox transformation.
  # Like the Classical Box Cox transformation but with a linear 
  # tail on the shrinking side, and a continuous derivative.
  #
  # Arguments:
  # y      : univariate data
  # lambda : transformation parameter for BC
  # chg    : change point: where the linear part begins
  #          when lambda < 1, or ends when lambda > 1.
  #          If NULL, uses the getChangepoint function
  # stdToo : also output poststandardized transformed data.
  #
  if(min(y) <= 0) stop("Data values should be strictly positive")
  yt = rep(NA,times = length(y))
  chgt <- NULL
  if(lambda == 1){ # is already linear
    yt = y-1
    chgt = 0
  }
  if (lambda > 1){
    if(is.null(chg)) { chg = getChangepointBCr(y, lambda = lambda) }
    indl = which(y < chg)
    if (length(indl) > 0) {
      yt[-indl] = (y[-indl]^lambda-1)/lambda
      chgt      = (chg^lambda-1)/lambda
      yt[indl]  = chgt + (y[indl]-chg)*chg^(lambda-1)
    } else {
      yt = (y^lambda-1)/lambda
    }
  } else if(lambda < 1){
    if(is.null(chg)) { chg = getChangepointBCr(y,  lambda = lambda) }
    indu = which(y > chg)    
    if (lambda == 0) {
      if (length(indu) > 0) {
        yt[-indu] = log(y[-indu])
        chgt      = log(chg)
        yt[indu]  = chgt + (y[indu]-chg)/chg
      } else {
        yt = log(y)
      }
    } else {
      if (length(indu) > 0) {
        yt[-indu] = (y[-indu]^lambda-1)/lambda
        chgt      = (chg^lambda-1)/lambda
        yt[indu]  = chgt + (y[indu]-chg)*chg^(lambda-1)
      } else {
        yt = (y^lambda-1)/lambda
      }
    }
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(yt, ncol = 1), 
                                        type = "hubhub")
      zt <- (yt - locScale$loc)/locScale$scale
    } else {
      zt <- yt
    }
  } else {
    zt = NULL
  }
  return(list(yt=yt, chg=chg, chgt=chgt, zt=zt))
}


YJr = function(y, lambda, chg = NULL, prec = 1e-10, stdToo = TRUE) {
  # Rectified Yeo-Johnson transformation.
  # Like the classical YJ transformation but with a linear tail 
  # on the shrinking side, and a continuous derivative.
  #
  # Arguments:
  # y      : univariate data. Assumed shifted to median zero.
  # lambda : transformation parameter for YJ.
  # chg    : change point: where the linear part begins
  #          when lambda < 1, or ends when lambda > 1.
  #          If NULL, uses Q3 for lambda > 1 and
  #          Q1 for lambda > 1.
  # stdToo : also output poststandardized transformed data.
  #
  quarts = localFivenum(y)[2:4]
  yt = rep(NA,times = length(y))
  if (lambda == 1){ # is already linear
    yt = y
    chgt = 0
  }
  indneg = which(y < 0)
  indpos = which(y >= 0)
  #
  if(lambda > 1){
    if(is.null(chg)) { chg = getChangepointYJr(y, lambda = lambda) }
    indl = which(y < chg) # to be linearized
    indb = which(y < 0 & y >= chg) # in between
    yt[indpos] = ((1 + y[indpos])^(lambda) - 1)/lambda
    #
    if(lambda == 2){
      yt[indb] = -log(1-y[indb])
      chgt     = -log(1-chg)
      yt[indl] = chgt + (y[indl]-chg) * YJprimex(chg,lambda)
    }
    else {
      yt[indb] = -((1-y[indb])^(2-lambda)-1)/(2-lambda)
      chgt     = -((1-chg)^(2-lambda)-1)/(2-lambda)
      yt[indl] = chgt + (y[indl]-chg) * YJprimex(chg,lambda)
    }
  }  
  if(lambda < 1){
    if(is.null(chg)) { chg = getChangepointYJr(y, lambda = lambda) }
    # if(chg <= 0) stop("chg should be positive")
    indu = which(y > chg) # to be linearized
    indb = which(y >= 0 & y <= chg) # in between
    yt[indneg] = -((1-y[indneg])^(2-lambda)-1)/(2-lambda)
    #
    if(lambda == 0){
      yt[indb] = log(1 + y[indb])
      chgt     = log(1 + chg)
      yt[indu] = chgt + (y[indu]-chg)*YJprimex(chg,lambda)
    }
    else {
      yt[indb] = ((1+y[indb])^(lambda)-1)/lambda
      chgt     = ((1+chg)^(lambda)-1)/lambda
      yt[indu] = chgt + (y[indu]-chg)*YJprimex(chg,lambda)
    }
  } 
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- estLocScale(matrix(yt, ncol = 1), 
                                        type = "hubhub")
      zt <- (yt - locScale$loc) / locScale$scale
    } else { 
      zt <- yt
    } 
  } else {
    zt <- NULL
  }
  return(list(yt=yt, chg=chg, chgt=chgt, zt=zt))
}

estML <- function(x, lambdarange = NULL, type = "BC", standardize = TRUE) {
  n <- length(x)
  if (is.null(lambdarange)) {
    lambdarange <- c(-4, 6)
  }
  if (type == "YJ") {
    locx <- 0
    scalex <- 1
    if (standardize) {
      locx <- mean(x)
      scalex <- sd(x)
      x <- (x - locx)/scalex
    }
    tempfunc <- function(lambdatemp) {
      xyj <- YJ(x, lambdatemp, stdToo = FALSE)$yt
      mu <- mean(xyj)
      sigma2 <- mean((xyj - mu)^2)
      result <- -0.5 * n * log(2 * pi) - 0.5 * n * log(sigma2) - 
        0.5/sigma2 * sum((xyj - mu)^2) + (lambdatemp - 
                                            1) * sum(sign(x) * log(1 + abs(x)))
      return(result)
    }
    opt.out <- optimize(tempfunc, range(lambdarange), maximum = TRUE)
    objective = tempfunc(opt.out$maximum)
    xt <- YJ(x, opt.out$maximum, stdToo = FALSE)$yt
    zt <- scale(xt)
  }
  else {
    locx <- 0
    scalex <- 1
    if (standardize) {
      scalex <- median(x)
      x <- x/scalex
    }
    tempfunc <- function(lambdatemp) {
      xyj <- BC(x, lambdatemp, stdToo = FALSE)$yt
      mu <- mean(xyj)
      sigma2 <- mean((xyj - mu)^2)
      result <- -n/2 * log(sigma2) + (lambdatemp - 1) * 
        sum(log(x))
      return(result)
    }
    opt.out <- optimize(tempfunc, range(lambdarange), maximum = TRUE)
    objective = tempfunc(opt.out$maximum)
    xt <- BC(x, opt.out$maximum, stdToo = FALSE)$yt
    zt <- scale(xt)
  }
  return(list(lambda = opt.out$maximum, lambdarange = lambdarange, 
              objective = objective, xt = xt, zt = zt, weights = rep(1, 
                                                                     n), locx = locx, scalex = scalex))
}

# ---------------------
#' print method for rewMLE class
#' @export
print.rewMLE  <- function(x, ...) {
  # Print methods do not need to use '...', 
  # but it's available if needed.
  
  # Print output
  cat("print class for 'rewMLE'\n\n")
  
  cat(str(x))
  
  invisible(x)
}

# ---------------------
#' plot method for rewMLE class
#' @export
plot.rewMLE <- function(x, trim = 0.95, ... , main = NULL){
  
  original <- x$X
  observed <- x$Y
  weights <- x$weights
  lambdahats <- x$lambdahats
  
  if (is.vector(observed)){ 
    original <- matrix(original, nrow = length(original), ncol = 1)
    observed <- matrix(observed, nrow = length(observed), ncol = 1)
    weights <- matrix(weights, nrow = length(weights), ncol = 1)
  }
  
  if (!identical(dim(observed),dim(weights))){
    stop("observed and weights need to have the same dimensions")
    
  } else if (!is.null(main) && length(main)!=dim(observed)[2]){
    stop("The names provided should be the same length as your data")
    
  } else {
    
    # Set the 'main' of plot
    if(is.null(main)){
      if(is.null(colnames(observed))){
        main <- seq_along(observed)
      } else {
        main <- colnames(observed)
      }
    }
    
    #par(mfrow = c(1, 2))
    layout(matrix(c(1,1,2,2,3,4,5,5), 2, 4, byrow = T))
    par(ps = 10)
    
    for(i in 1:ncol(observed)){
      
      x_orig <- original[,i]
      x_i <- observed[,i]
      y_i <- seq_along(x_i)
      w_i <- weights[,i]
      
      name <- main[i]  
      
      #par(mfrow = c(1, 1), ps = 10)
      
      #######################
      # Regular histogramm
      # ::::::::::::::::::::
      
      breaks_i <- seq(min(x_orig, na.rm = T),
                      max(x_orig, na.rm = T),
                      by = diff(range(x_orig, na.rm = T))/50)
      
      # Plot histogram with weights == 1
      hist(x_orig[w_i==1], breaks = breaks_i, col = rgb(1, 0, 0, 0.5), 
           main = paste("Original data (+ suspected outliers): ", name), 
           xlab = "Value", 
           xlim = range(x_orig, na.rm = T)#, 
           #ylim = c(0, max(hist(x_orig[w_i==1], plot = F)$counts))
      )
      
      # Overlay histogram with weights == 0
      hist(x_orig[w_i==0], 
           breaks = breaks_i, 
           col = rgb(0, 0, 1, 0.5), 
           add = TRUE)
      
      #######################
      # Highlighting outliers
      # TRUE: Non-outlier 
      # FALSE: Outlier
      # -----------------
      
      
      # Plotting the datapoints and highlighting the detected outliers
      plot(x_i, y_i, type='n', 
           main = paste("Transformed data with outliers: ", name), 
           xlab = "Value", ylab = "Rank") 
      
      # Adding points based on condition
      points(x_i[w_i==1], y_i[w_i==1], col='black') # Points where TRUE
      points(x_i[w_i!=1], y_i[w_i!=1], col='red', pch=8) # Points where FALSE
      
      abline(v = c(min(x_i[w_i==1], na.rm = T), max(x_i[w_i==1], na.rm = T)), 
             col = "black", lty = 3 )
      
      # Adding a legend
      legend("topright", legend=c("Outlier", "Non-outlier"), 
             col=c("red", "black"), pch=c(8, 1))
      
      
      
      ##################
      # Plotting the CDF
      # ----------------
      
      #par(mfrow = c(1, 2), ps = 10)
      
      # Generate QQ plot data for x
      qq <- qqnorm(x_i, plot = FALSE)  
      
      # For trim equals NULL
      # --> estimate which portion would minimally cover it
      if (is.null(trim)){
        trim <- estTrim(data = qq$x, weights = w_i)
      }
      
      # Trimmed R-squared value for labeling
      R_trim <- mseTrim(x_i, trim = trim)
      Qts <- qnorm(c((1-trim)/2, 1-(1-trim)/2))
      
      # Setup the plot
      plot(qq$x, qq$y, type='n', 
           main=bquote("QQ: " ~ R^2 ~ " of "
                       ~ .(round(trim,3)*100) ~ "%: " ~ .(round(R_trim,2))), 
           xlab="Theoretical Quantiles", ylab="Sample Quantiles")
      
      points(qq$x[w_i == 1], qq$y[w_i == 1], col='black')  # Normal points
      
      # Check if there are any outliers: If yes, color them red
      if(any(w_i != 1)) {  
        points(qq$x[w_i != 1], qq$y[w_i != 1], col='red', pch=8)
      }
      
      qqline(x_i, col = "black")
      abline(v = Qts, col = "black", lty = 3)
      
      
      # Residual Plot
      x_i <- x_i[w_i == 1]
      
      n <- length(x_i)
      probs <- (rank(x_i) - 0.5) / n
      # Expected values if x from a normal distr.
      expected <- qnorm(probs)  
      
      # Calculate residuals as difference from observed to expected
      residuals <- x_i - expected
      plot(x_i, residuals, xlab = "Observed Value", 
           ylab = "Residual", 
           main = "Normal residuals plot")
      
      abline(h = 0, col = "black")
      abline(v = Qts, col = "black", lty = 3)
      
      #######################
      # Transformed histogramm
      # ::::::::::::::::::::
      
      breaks_j <- seq(min(x_i, na.rm = T),
                      max(x_i, na.rm = T),
                      by = diff(range(x_i, na.rm = T))/50)
      
      # Plot histogram with w_i == 1
      hist(x_i[w_i==1], breaks = breaks_j, col = rgb(1, 0, 0, 0.5), 
           main = bquote("Transformed data: " ~ R^2 ~ "=" ~ .(round(R_trim, 3)) ~ "," 
                         ~ hat(lambda) ~ "=" ~ .(round(lambdahats[i], 2))),
           xlab = "Value")
      #xlim = range(x_i, na.rm = T))
      
      abline(v = Qts, col = "black", lty = 2)
      
      # # Overlay histogram with w_i == 0
      # hist(x_i[w_i!=1], breaks = breaks_j, col = rgb(0, 0, 1, 0.5), 
      #      add = TRUE)
      
      
      
    }
  }
}


#' Estimate the central maximal Interval (eg. 95%) that
#' @export
estTrim <- function(data, weights){
  valid <- !is.na(data) & weights == 1
  data <- data[valid]
  
  z_val_min <- min(data)
  z_val_max <- max(data)
  
  # Smallest percentile that still fits (on the left)
  perc_min <- pnorm(z_val_min)
  
  # Highest percentile that still fits (int the rights)
  perc_max <- 1 - pnorm(z_val_max)
  
  # Which percentile is closer to the mean? 
  # Meaning we are looking for the bigger value
  perc_closer_to_mean <- ifelse(perc_min > perc_max, 
                                perc_min, perc_max)
  
  # Return the 0.xx range (xx% Interval)
  return(1-2*perc_closer_to_mean)
  
}


#' Mean Square Error (MSE) of a trimmed distribution
#' 
#' This function assesses how closely a trimmed version 
#' of the input data follows a normal distribution by 
#' comparing the observed values against expected values 
#' derived from a normal distribution fitted to the data. 
#' It returns the R^{2} value as a measure of this fit. 
#' @export
mseTrim <- function(data, trim = 0.95){
  
  if(!is.numeric(trim) || trim < 0 || trim > 1) {
    stop("trim must be a numeric value between 0 and 1")
  }
  
  obs <- sort(data, na.last = NA)
  p <- (1:length(obs) - 0.5) / length(obs)
  exp <- qnorm(p = p, mean = mean(obs), sd = sd(obs))
  
  # trimming
  lower <- which.min(p < (1-trim)/2)
  higher <- which.max(p > 1-(1-trim)/2) - 1
  
  # Checking for valid indices
  lower <- max(1, lower)
  higher <- min(length(obs), higher)
  
  obs <- obs[lower:higher]
  exp <- exp[lower:higher]
  
  # Calculate Mean Standard Error
  mse <- round(mean((obs - exp)^2, na.rm = T), digits= 3)
  # Calculate correlation and return the R^2 value
  corr <- round(cor(obs, exp)^2, digits = 5)
  return(corr)
  
}

