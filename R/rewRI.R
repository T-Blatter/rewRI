############################################################
## RI functions: estimating RI by re-weighted MLE (BC/YJ) ##
############################################################

#####################
#' Estimating the RIs from a previously established re-weighted MLE object
#' 
#' This function takes as an input either an object of class 
#' "rewMLE" or "rewMLE.boot". For estimation of RIs diretly from
#' clinical data, refer to the `estRI` function, which is extended
#' wrapper function of this function.
#' @export
rewRI <- function(rewMLE.out, RI = 0.95){
  
  if (!inherits(rewMLE.out, c("rewMLE", "rewMLE.boot"))){
    stop("rewMLE.out should either be of class \"rewMLE\" or \"rewMLE.boot\".")
  }
  
  # Single estimate
  if(inherits(rewMLE.out, "rewMLE")){
    # Data
    X <-          rewMLE.out$X
    # Data: tf. & std.
    Y <-          rewMLE.out$Y
    # Data: weights
    w <-          rewMLE.out$weights
    
    # estimated Parameters
    lambdahats <- rewMLE.out$lambdahats
    muhat <-      rewMLE.out$muhat
    sigmahat <-   rewMLE.out$sigmahat
    
    # Quantiles
    qthatSim  <-          matrix(NA, ncol=ncol(Y), nrow = 2)
    qthatEst  <-          matrix(NA, ncol=ncol(Y), nrow = 2)
    colnames(qthatEst) <- colnames(qthatSim) <- colnames(Y)
    rownames(qthatEst) <- rownames(qthatSim) <- c("LL", "UL")
    
    for (j in 1:ncol(Y)){
      
      y <- Y[,j][as.logical(w[,j])]
      
      # --- 
      # Regular Quantiles
      qthatSim[,j] <-   qnorm(p = c((1-RI)/2, 1-(1-RI)/2),
                              mean = muhat[j], sd = sigmahat[j] )
      
      # ---
      # TruncNormal Quantiles
      l.limit <- min(y, na.rm = T)
      u.limit <- max(y, na.rm = T)
      
      # Calculate the proportions of the 
      # original distribution outside the truncation limits
      cdfs <-  pnorm(c(l.limit,u.limit), mean = muhat[j], 
                     sd = sigmahat[j])
      
      # Adjust your quantile levels for the truncation
      # This adjusts the quantile levels to the 
      # proportion of the distribution that remains
      adjusted_qt <- cdfs[1] + (cdfs[2] - cdfs[1]) * c((1-RI)/2, 1-((1-RI)/2))
      
      # Use the adjusted quantiles with the inverse 
      # CDF of the normal distribution
      qthatEst[,j] <- qnorm(adjusted_qt, mean = muhat[j], sd = sigmahat[j])
      
    }
    
    qtSim  <- transfoback_rewMLE(Ynew = qthatSim, rewMLE.out = rewMLE.out)
    qthat  <- transfoback_rewMLE(Ynew = qthatEst, rewMLE.out = rewMLE.out)
    
    # Print the quantiles
    for (j in 1:ncol(Y)){
      
      cat("     Norm Quantiles:\t", qtSim[1,j], "\t", qtSim[2,j], "\n")
      cat("truncNorm Quantiles:\t", qthat[1,j], "\t", qthat[2,j], "\n\n")
      
    }
    
    result    <- c(list("RI" = RI,
                        "boot" = F,
                        # Qts in BC/YJ space
                        "qthatEst" = qthatSim, #qthatEst
                        # Qts in non-tf space
                        "qthat" = qtSim), #qthat
                   rewMLE.out)
    
  } 
  
  # Bootstrapping estimates
  if(inherits(rewMLE.out, "rewMLE.boot")){
    
    # How many times?
    nBoot <-      rewMLE.out$nBoot
    
    # Original Data
    X <-            rewMLE.out$X
    dorig <-        rewMLE.out$dorig
    standardize <-  rewMLE.out$standardize
    ttypes <-       rewMLE.out$ttypes
    
    # estimated Parameters
    lambdahats <- rewMLE.out$lambdahats
    muhat <-      rewMLE.out$muhat
    sigmahat <-   rewMLE.out$sigmahat
    locx <-         rewMLE.out$locx
    scalex <-       rewMLE.out$scalex
    
    # Weights
    w <-          rewMLE.out$weights
    
    # Quantiles
    qtBoot_LL <-        matrix(NA, ncol=ncol(X), nrow = nBoot)
    qtBoot_UL <-        matrix(NA, ncol=ncol(X), nrow = nBoot)
    
    colnames(qtBoot_LL) <- colnames(qtBoot_UL) <- colnames(X)
    
    # Iterate over n bootstraps
    for (i in 1:nrow(lambdahats)){
      
      tf_list <- list(dorig=dorig,
                      standardize=standardize,
                      locx=locx[i,],
                      scalex=scalex[i,],
                      muhat=muhat[i,],
                      lambdahats=lambdahats[i,],
                      sigmahat=sigmahat[i,],
                      ttypes=ttypes[i,])
      
      Y <- transfo_rewMLE(X, rewMLE.out = tf_list)
      
      qtSim <- qts <- matrix(NA, ncol=ncol(Y), nrow = 2)
      
      
      for (j in 1:ncol(Y)){
        
        y <- Y[,j][as.logical(w[,j])]
        
        # -- 
        # Normal Quantiles
        qtSim[,j] <-    qnorm(p = c((1-RI)/2, 1-(1-RI)/2),
                              mean = muhat[i,j], sd = sigmahat[i,j] )
        # ---
        # TruncNormal Quantiles
        l.limit <- min(y, na.rm = T)
        u.limit <- max(y, na.rm = T)
        
        # Calculate the proportions of the 
        # original distribution outside the truncation limits
        cdfs <- pnorm(c(l.limit,u.limit), mean = muhat[i,j], sd = sigmahat[i,j])
        
        # Adjust your quantile levels for the truncation
        # This adjusts the quantile levels to the 
        # proportion of the distribution that remains
        adjusted_qt <- cdfs[1] + (cdfs[2] - cdfs[1]) * c((1-RI)/2, 1-((1-RI)/2))
        
        # Use the adjusted quantiles with the inverse 
        # CDF of the normal distribution
        qts[,j] <- qnorm(adjusted_qt, mean = muhat[i,j], sd = sigmahat[i,j])
        
        
      }
      
      qtSim_back  <-   transfoback_rewMLE(qtSim, rewMLE.out = tf_list)
      qts_back    <-   transfoback_rewMLE(qts,   rewMLE.out = tf_list)
      
      qtBoot_LL[i,] <- qtSim_back[1, ] #qts_back[1, ]
      qtBoot_UL[i,] <- qtSim_back[2, ] #qts_back[2, ]
      
    }
    
    # Quantiles: Calculate the Bootstrapped CIs
    qthat <-           matrix(NA, ncol=ncol(Y), nrow = 6)
    rownames(qthat) <- c("LL", "UL", "LL.low", "LL.high", "UL.low", "UL.high")
    colnames(qthat) <- colnames(Y)
    
    
    qtBoot_LL[is.infinite(qtBoot_LL) | is.nan(qtBoot_LL)] <- NA
    qtBoot_UL[is.infinite(qtBoot_UL) | is.nan(qtBoot_UL)] <- NA
    
    qtMeans_LL <- colMeans(qtBoot_LL, na.rm = T)
    qtMeans_UL <- colMeans(qtBoot_UL, na.rm = T)
    
    CIn <- matrix(NA, ncol=ncol(Y), nrow = 4)
    
    # Bias Corrected and Accelerated CIs 
    # BCa <- matrix(NA, ncol=ncol(Y), nrow = 4)
    # for (k in 1:ncol(Y)){
    #   BCa[1:2,k] <- estBCa(qtBoot_LL, orig_stat = qtMeans_LL[k], CI=0.9)
    #   BCa[3:4,k] <- estBCa(qtBoot_UL, orig_stat = qtMeans_UL[k], CI=0.9)
    # }
    
    CIn[1:2,] <-  apply(X = qtBoot_LL, 
                        MARGIN = 2, 
                        FUN = quantile, p=c(0.05,0.95), na.rm=T)
    
    CIn[3:4,] <-  apply(X = qtBoot_UL, 
                        MARGIN = 2, 
                        FUN = quantile, p=c(0.05,0.95), na.rm=T)
    
    # Assign calculated values
    qthat[1, ]    <- qtMeans_LL
    qthat[2, ]    <- qtMeans_UL
    qthat[3:6, ]  <- CIn #BCa
    
    qthatEst <- transfo_rewMLE(qthat,
                               rewMLE.out = list(dorig=dorig,
                                                 standardize=standardize,
                                                 locx=colMeans(locx),
                                                 scalex=colMeans(scalex),
                                                 muhat=colMeans(muhat),
                                                 lambdahats=colMeans(lambdahats),
                                                 sigmahat=colMeans(sigmahat),
                                                 ttypes=ttypes[1,]))
    
    result        <- c(list("RI" = RI,
                            "boot" = T,
                            # Qts in BC/YJ space
                            "qthatEst" = qthatEst,
                            # Qts in non-tf space
                            "qtBoot_LL" = qtBoot_LL,
                            "qtBoot_UL" = qtBoot_UL,
                            "qthat" = qthat),
                       rewMLE.out)
    
  }
  
  class(result) <- "rewRI"
  return(result)
  
  
  
}
# -----------------------
#' Printing method for the `rewRI` function
print.rewRI <- function(x){
  
  cat(str(x))
  
  invisible(x)
}


# -----------------------
#' Plotting method for the `rewRI` function
plot.rewRI <- function(x){

  RI    <- x$RI
  X     <- x$X
  Y     <- x$Y
  w     <- x$weights

  # Model Params
  lambdahats <- x$lambdahats
  muhat      <- x$muhat
  sigmahat   <- x$sigmahat

  # Quantiles
  qthat     <- x$qthat
  qthatEst  <- x$qthatEst

  # Density estimate from the data
  N  <-   50
  Xint  <-      matrix(NA, ncol=ncol(Y), nrow = N)
  freqX   <-    matrix(NA, ncol=ncol(Y), nrow = N)
  densityX <-   matrix(NA, ncol=ncol(Y), nrow = N)
  colnames(Xint) <- colnames(freqX) <- colnames(Y)
  colnames(densityX) <- colnames(Y)

  # X Scale estimate from the data (BC/YJ Space)
  Yint     <-   matrix(NA, ncol=ncol(Y), nrow = N)
  freqY   <-    matrix(NA, ncol=ncol(Y), nrow = N)
  densityY   <- matrix(NA, ncol=ncol(Y), nrow = N)
  colnames(Yint) <- colnames(freqY) <- colnames(Y)
  colnames(densityY) <- colnames(Y)

  # Density estimate from the rewMLE output (BC/YJ Space)
  Yest     <-   matrix(NA, ncol=ncol(Y), nrow = N)
  freqYest   <-  matrix(NA, ncol=ncol(Y), nrow = N)
  densityYest <- matrix(NA, ncol=ncol(Y), nrow = N)
  colnames(Yest) <- colnames(freqYest) <- colnames(Y)
  colnames(densityYest) <- colnames(Y)

  # Density estimates back-tf from tf space
  Xorg <-         matrix(NA, ncol=ncol(Y), nrow = N)
  freqXorg   <-   matrix(NA, ncol=ncol(Y), nrow = N)
  densityXorg <-  matrix(NA, ncol=ncol(Y), nrow = N)
  colnames(Xorg) <- colnames(freqXorg) <- colnames(Y)
  colnames(densityXorg) <- colnames(Y)

  # Transform X to BC/YJ Space: X.y
  X.y    <-  transfo_rewMLE(Xnew = X, rewMLE.out = x)
  # Transform Y in BC/YJ Space: Y.x
  Y.x    <-   transfoback_rewMLE(Ynew = Y[w], rewMLE.out = x)

  par(mfrow = c(2, 2), cex.main = 0.8)

  for(i in 1:ncol(X)){


    # Regular Space: Data
    Xmin <- min(X[,i], na.rm = T)
    Xmax <- max(X[,i], na.rm = T)
    Xbreaks   <-    seq(Xmin, Xmax, length.out = N+1)
    Xhist     <-    hist(x = X[,i], breaks = Xbreaks, plot = F)
    Xint[,i]  <-    Xhist$mids
    freqX[,i]   <-  Xhist$counts
    densityX[,i] <- Xhist$density

    # BC/YJ Space: Data
    Ymin <- min(X.y[,i], na.rm = T)
    Ymax <- max(X.y[,i], na.rm = T)
    Ybreaks   <-    seq(Ymin, Ymax, length.out = N+1)
    Yhist     <-    hist(x = X.y[,i], breaks = Ybreaks, plot = F)
    Yint[,i]  <-    Yhist$mids
    freqY[,i]   <-  Yhist$counts
    densityY[,i] <- Yhist$density

    # BC/YJ Space: rewMLE estimates
    Yweight   <-        Y[,i][w[,i]==1]
    Yest[,i]  <-        Yint[,i]
    freqYest  <-        dnorm(x = Yest[,i], mean = muhat[i],
                              sd = sigmahat[i]) * length(Yweight)
    densityYest[,i] <-  density(Yweight, bw = "nrd0", n = N, na.rm=T)$y

    # Regular Space: rewMLE estimates
    Xhistorg  <-        hist(Y.x[,i], breaks = Xbreaks, plot = F)
    Xorg[,i] <-         Xhistorg$mids
    freqXorg[,i]   <-   Xhistorg$counts
    densityXorg[,i]  <- Xhistorg$density

    plot(Xhist, main = "Normal Space", fill = "gray30")
    plot(Xhistorg, add=T, fill= "gray44")

    plot(Yhist, main = "Normal Space", fill = "gray30")
    plot(Xhistorg, add=T, fill= "gray44")

    # # Transformed SPACE
    # plot(Yint[,i], densityY[,i],
    #      main = paste0("Transformed space: ", colnames(Y)[i], " with lambda ",
    #                    round(lambdahats[i],2)),
    #      xlab = "Value", ylab = "Density", pch=19, cex=0.2, col = "darkgreen")
    #
    # lines(Yint[,i], densityY[,i], lty = 1, lwd = 1, col = "darkgreen")
    #
    # points(Xorg[,i], densityXorg[,i], pch=19, cex=0.2, col = "lightblue")
    # lines(Xorg[,i], densityXorg[,i], lty = 1, lwd = 1, col = "lightblue")
    #
    # abline(v=qthat[,i], col= "green4", lty = 2)
    #
    # # Original SPACE
    # plot(XEst[,i], densityY[,i],
    #      main = paste0("Original space: ", colnames(Y)[i]),
    #      xlab = "Value", ylab = "Density", pch=19, cex=0.2, col = "darkgreen")
    # lines(XEst[,i], densityY[,i], lty = 1, lwd = 0.5, col = "darkgreen")
    #
    # points(XorgEst[,i], densityXorg[,i], pch=19, cex=0.2, col ="lightblue")
    # lines(XorgEst[,i], densityXorg[,i], lty = 1, lwd = 0.5, col ="lightblue")
    #
    # abline(v=qthatEst[,i], col= "green4", lty = 2)

  }

}
