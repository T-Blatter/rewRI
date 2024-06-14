############################################################
## RI functions: estimating RI by re-weighted MLE (BC/YJ) ##
############################################################

#####################
#' Estimating the RLs
#' 
#' This is a wrapper function, that estimates RIs either directly from the 
#' data or from a previously established re-weighted RI.
#' @export
estRI <- function(X, RI=0.95, nBoot=NULL, ...){
  
  # Input: X = Data (no class)
  if (!inherits(X, c("rewMLE", "rewMLE.boot", "rewRI"))) {
    
    x <- prepData(X)
    
    # Single Estimate
    if (is.null(nBoot)){
      rewMLE.out    <-     rewMLE(X = x$X, ...)
      
    # Bootstrapping with nBoot reps
    } else {
      rewMLE.out   <-      rewMLE.boot(X = x$X, nBoot=nBoot, ...)
      
    }
  
  # Input: X = Class "rewMLE" / "rewMLE.boot"
  } else if (!inherits(X, "rewRI")) {
    rewMLE.out <-          X
    
    
  # Input: X = Class "rewRI"
  } else {
    result <- list(X,
                   Data=x$X)
    
    class(result) <- "estRI"
    return(result)
  }
  
  
  # Estimate RIs given the rewMLE.out
  result        <- rewRI(rewMLE.out, RI = RI)
  class(result) <- "estRI"
  return(result)
  
}

#' Prepare data for re-weighted MLE
prepData <- function(Data){
  
  # Check if vector or matrix
  if (is.vector(Data)) {
    X <- matrix(Data, ncol = 1)
  } 
  
  X <- as.matrix(Data)
  W <- as.matrix(Data)
  W[T] <- 1
  n <- nrow(X)
  d <- ncol(X)
  
  skew <- rep(NA, d)
  
  for (i in seq_len(d)) {
    
    #########
    # calculate Skewness
    skew[i] <- estSkew(X[,i])
    
    #########  
    # # median +/- two median absolute deviations (MAD)
    # two_mad_range <- c(median(X[,i], na.rm=T) - 2*mad(x = X[,i], na.rm=T),
    #                    median(X[,i], na.rm=T) + 2*mad(x = X[,i], na.rm=T))
    # 
    # W[X[,i] < two_mad_range[1] | X[,i] > two_mad_range[2], i] <- 0
    
    # # mean +/- three sds
    # three_sigma_range <- c(mean(X[,i], na.rm=T) - 2*sd(x = X[,i], na.rm=T),
    #                        mean(X[,i], na.rm=T) + 2*sd(x = X[,i], na.rm=T))
    # W[X[,i] < three_sigma_range[1] | X[,i] > three_sigma_range[2], i] <- 0
    
  }
  
  return( list("X" = X,
               "W" = W,
               "skew" = skew))
  
  
}

#' Estimate Skewness from the Data
estSkew <- function (x, na.rm=TRUE) {
  
  if (is.matrix(x)) { apply(x, 2, estSkew, na.rm=na.rm)  }
  
  else if (is.vector(x)) {
    if (na.rm) { x <- x[!is.na(x)] }
    
    n <- length(x)
    return( (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2) )
    
  } else { estSkew(as.vector(x), na.rm=na.rm) }
  
}

#' Reference Interval from a Normal Distribution
#' 
#' Calculate the Reference Interval (+ 90% Confidence Interval) from a normal distribution 
#' given the population mean (mu) and standard deviation (sd, sigma)
#' 
#' Confidence intervals of sample error
#' se = sqrt(((sd^2)/n) + (((z.score.RI^2)*(sd^2))/(2*n)))
#' 
#' https://en.wikipedia.org/wiki/Reference_range#Confidence_interval_of_limit
tnormRI <- function(n = 1000, mean = 0, sd = 1, RI = 0.95, CI = 0.9){
  
  RIbounds <- c((1 - RI) / 2, (1 + RI) / 2)
  CIbounds <- c((1 - CI) / 2, (1 + CI) / 2)
  
  RIs <- qnorm(p = RIbounds, mean = mean, sd = sd)
  
  z.score.RI <- qnorm(1 - ((1 - RI) / 2))
  z.score.CI <- qnorm(1 - ((1 - CI) / 2))
  
  se <- sqrt((sd^2/n) + (z.score.RI^2*(sd^2))/(2*n))
  
  return( c(RIs[1], RIs[2],
            RIs[1] - z.score.CI * se, 
            RIs[1] + z.score.CI * se,
            RIs[2] - z.score.CI * se,
            RIs[2] + z.score.CI * se ) )
}


#' print function of estRI class
#' @export
print.estRI <- function(x){
  
  cat("Estimated Reference Interval\n")
  if(x$boot){ cat(paste0(x$nBoot , " bootstrap iterations:\n"))}
  cat(">--------------------------<\n\n")
  
  RI <- x$RI
  qts <- c((1-RI)/2,1-((1-RI)/2))
  
  cat(paste0(RI*100, "% RI\t", paste0(c((1-RI)/2,1-((1-RI)/2)), collapse = "\t"),
             ifelse(x$boot, yes="\t90%CIs\n\n" , no="\n\n")))
  
  estRI <- round(x$qthat, 3)
  estRInames <- colnames(estRI)
  estQtnames <- rownames(estRI)
  
  cat("\t", paste(estQtnames, collapse = "\t"),"\n")
  
  for(i in 1:ncol(estRI)) {
    cat(estRInames[i], "\t", paste(estRI[,i], collapse = "\t"), "\n", sep="")
  }
  
  cat(">--------------------------<\n\n")
  
  # Custom return
  invisible(x)
}

#' plot function of estRI class
#' @export
plot.estRI <- function(x, xname = NULL){
  
  data  <- x$X
  RI    <- x$RI
  
  # Density from data
  Xorg  <- x$Xorg
  Yorg  <- x$Yorg
  
  # Estimates in transported space
  X         <- x$X
  Y         <- x$Y
  densityY  <- x$densityY
  qthat     <- x$qthat
  
  # Estimates back in original space
  XorgEst   <- x$XorgEst
  XEst      <- x$XEst
  qthatEst  <- x$qthatEst
  
  # Lambda estimates
  lambdahats <- x$lambdahats
  
  names <- colnames(data)
  
  par(mfrow = c(1, 1))
  
  for (j in 1:ncol(data)){
    
    d <- data[,j]
    d <- d[!is.na(d)]
    
    # Histogram with Original Data
    hist(d, breaks = 100, probability = TRUE, col = "white", border = "grey",
         main = paste0(names[j], " data: estimated ", 
                       RI*100, "% RI"),
         ylab = "Density", xlab = "Value")
    
    # Density Overlays
    D_dens <- density(d, bw = "nrd0", na.rm=T)
    
    # Estimated density
    X_seq <- XEst[,j]
    Y_dens <-  Yorg[,j]
    Y_densEst <- densityY[,j]
    
    # Scale factor
    scale_fact <- (max(D_dens$y) / max(Y_dens))
    
    Y_dens <- Y_dens * scale_fact
    Y_densEst <- Y_densEst * scale_fact
    
    # Actual density and estimated density curve
    #lines(XorgEst[,j], Y_dens, col = "lightblue", lwd = 2)
    lines(D_dens, col ="lightblue4", lwd = 2)
    
    points(XorgEst[,j], Y_dens, pch=19, cex=0.2, col ="lightblue")
    lines(XorgEst[,j], Y_dens, lty = 1, lwd = 0.5, col ="lightblue")
    
    lines(X_seq, Y_densEst, col = "lightgreen", lwd = 2, lty = 1) 
    
    # Add Vertical Lines for Reference Percentiles
    RIs <- x$qthat[,j]
    abline(v = RIs, col = "green4", lwd = 2, lty = 2)
    
    legend("topright", legend = c("Data Density", "Estimated Density", "Reference Interval"), 
           col = c("lightblue4", "darkgreen", "green4"), lwd = 2, lty = c(1, 1, 2))
    
  }
  
}
