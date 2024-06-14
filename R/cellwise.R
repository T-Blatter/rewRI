##############
# Functions from the cellWise package
# cellWise: Analyzing Data with Cellwise Outliers
# Version: 2.5.2
# License: GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]
# https://cran.r-project.org/web/packages/cellWise/index.html

localFivenum = function(x, na.rm = TRUE) {
  # Our version of stats::fivenum where now the quartiles
  # are order statistics, with symmetric ranks.
  # Assumes x is _sorted_ already.
  xna <- is.na(x)
  if (any(xna)) {
    if (na.rm) 
      x <- x[!xna]
    else return(rep.int(NA, 5))
  }
  n <- length(x)
  if (n == 0) 
    rep.int(NA, 5)
  else {
    med <- if ((n %% 2) == 0) {
      (x[(n/2)] + x[(n/2) + 1]) / 2
    } else {
      x[(n - 1)/2 + 1]
    }
    result <- c(x[1], x[ceiling(n / 4)], med,
                x[n - ceiling(n / 4) + 1], x[n])
    result
  }
}


rhoBiweight <- function(x, c = 0.5) {
  xbw <- c * (1 - (1 - (x/c)^2)^3)
  xbw[which(abs(x) > c)] <- c
  return(xbw)
}


robnormality <- function(y, b = 0.5) {
  # Assess normality of a dataset after robustly standardizing.
  # Outliers only have a limited effect on this criterion.
  # y is assumed to be _sorted_.
  #
  # args: 
  #   y : vector of observations
  #   b : tuning parameter of biweight rho
  # 
  sy       <- y[!is.na(y)]
  n        <- length(sy)
  locScale <- estLocScale(matrix(sy, ncol = 1), 
                          type = "hubhub")
  sy       <- (sy - locScale$loc) / locScale$scale
  # Theoretical quantiles and differences:
  hprobs     <- ((seq_len(n)) - 1/3)/(n + 1/3)
  theoQs     <- qnorm(hprobs)
  diffs      <- sy - theoQs
  # apply Biweight rho
  rhodiffs <- rhoBiweight(diffs, b)
  crit     <- mean(abs(rhodiffs))
  return(crit)
}


estLocScale <- function(X, 
                        type = "wrap", precScale = 1e-12,
                        center = TRUE,
                        alpha = 0.5,
                        nLocScale = 25000,
                        silent = FALSE) {
  # Estimates location and scale for every column in X
  #
  # args:
  #   X: data matrix
  #   nLocScale: if < n, then rows are sampled tocalculate loc/scale
  #   type: "wrap" or "mcd", the location/scale estimators used
  #   precScale: precision scale used throughout the algorithm
  #   alpha = value of alpha for the unimcd, h = ceil(alpha * n)
  # Returns: 
  #   loc: the locations of the columns in X
  #   scale: the scales of the columns in X
  
  # Check inputs
  
  # The random seed is retained when leaving the function
  if (exists(".Random.seed",envir = .GlobalEnv,inherits = FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv, 
                     inherits = FALSE)
    on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
  }
  set.seed(0)
  
  if (!is.data.frame(X) & !is.matrix(X) & !is.vector(X)) {
    stop("The input data must be a vector, matrix or a data frame")
  }
  type <- match(type, c("1stepM","hubhub","wrap","mcd", "rawmcd", "wrapmedmad")) - 1
  if (is.na(type)) {
    stop(paste("Invalid \"type\" argument. Should be \"wrap\", \"mcd\" or \"1stepM\""))
  }
  
  if (is.na(alpha)) {
    alpha <- 0.5
  }
  
  if (nLocScale == 0) {
    nLocScale = dim(X)[1]
  }
  
  # Estimate location/scale
  res <- tryCatch( .Call('_cellWise_estLocScale_cpp', as.matrix(X), nLocScale,
                         type,  precScale,
                         center, alpha, PACKAGE = 'cellWise'),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  zeroscales <- which(res$scale <= precScale)
  if (!silent) {
    if ( length(zeroscales) > 0) {
      warning(paste(length(zeroscales)," out of ", dim(X)[2], " variables have an estimated scale <= 
\"precScale\" = ", precScale, "."))
    }
  }
  
  loc.out <- drop(res$loc)
  scale.out <- drop(res$scale)
  names(loc.out) <- names(scale.out) <- colnames(X)
  return(list(loc = loc.out, scale = scale.out))
}


getCritval <- function(x, lambda, tfunc = YJr, chg2=NULL,  quant=0.99) {
  # Calculates the value of the robust fitting criterion.
  # Arguments:
  # x      : the univariate data
  # lambda : the value of lambda
  # tfunc  : either BC, BCr, YJ, YJr
  # chg2   : 2 changepoints: one for lambda >= 1 and one if lambda < 1
  #          so chg[1] < chg[2]
  #          If NULL, uses Q1 for lambda > 1 and Q3 for lambda > 1
  #
  if (is.null(chg2)) { chg = NULL } else {
    chg <- if (lambda >= 1) {chg2[1]} else {chg2[2]}
  }
  robnormality(tfunc(x,lambda,chg,stdToo = FALSE)$yt, b = 0.5)
}


getChangepointYJr <- function(y, lambda, fac = 1.5, eps = 1e-5) {
  # Function to caculate the "changepoint" for the rectified
  # Yeo-Johnson transformation
  # 
  quarts = localFivenum(y)[2:4]
  if (lambda < 1){
    chg = YJ(quarts[3], lambda, stdToo = FALSE)$yt * 1.5
  } else if (lambda > 1){
    chg = YJ(quarts[1], lambda, stdToo = FALSE)$yt * 1.5
  }
  if (lambda < 0) {
    chg <- min(chg, abs(1 / lambda) - eps)
  } else if (lambda > 2) {
    chg <- max(chg, (1 / (2 - lambda)) + eps)
  }
  chg <- iYJ(chg, lambda, stdToo = FALSE)$yt
  chg <- min(max(chg, y[1]), y[length(y)])
  return(chg)
}


getChangepointBCr <- function(y, lambda, fac = 1.5, eps = 1e-5) {
  # Function to caculate the "changepoint" for the rectified
  # Box Cox transformation
  # 
  quarts = localFivenum(y)[2:4]
  if (lambda < 1){
    chg = BC(quarts[3], lambda, stdToo = FALSE)$yt * 1.5
  } else if (lambda > 1){
    chg = BC(quarts[1], lambda, stdToo = FALSE)$yt * 1.5
  }
  #correct for image of BC
  if (lambda < 0) {
    chg <- min(chg, abs(1 / lambda) - eps)
  } else  if (lambda > 0) {
    chg <- max(chg, -(1 / lambda) + eps)
  }
  chg <- iBC(chg, lambda, stdToo = FALSE)$yt
  chg <- min(max(chg, y[1]), y[length(y)])
  return(chg)
}


getTtype <- function(x, type, j=NULL) {
  # Gets the transformation type (BC, YJ, or bestObj) and checks 
  # whether the domain of the variable is compatible with that type.
  # args:
  #   x    : vector with the values of the variable
  #   type : one of "BC", "YJ", "bestObj"
  #   j    : number of the variable
  # 
  if (!type %in% c("BC", "YJ", "bestObj")) {
    stop(" type should be 'BC' or 'YJ' or 'bestObj'.")
  }
  minx <- min(x)
  if (type == "BC") {
    if (minx <= 0) {
      stop(paste("The minimum of variable ", j,
                 " is not strictly positive.\n",
                 "Please choose a type other than BC.",sep = ""))
    } else {
      ttype <- "BC"
    }
  } else if (type == "YJ") {
    ttype <- "YJ"
  } else if (type == "bestObj") {
    if (minx <= 0) {
      ttype <- "YJ"
    } else {
      ttype <- "bestObj"
    }
  }
  return(ttype)
}


