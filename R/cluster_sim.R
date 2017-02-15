
#' Calculate expected coverage of several alternative intervals for
#'  marginal coefficients in linear regression with clustered coefficients.
#' @description This function can be used to generate simulations in section 3.1
#' @param beta List of regression coefficients. 
#' This should be a least with K elements where K is the number of blocks
#' @param which.sample Indices of individuals who are sampled
#' @param Sigma List of correaltion matrices, one for each block
#' @return A list with elements:
#' \describe{
#' \item{\code{COVERAGE}}{A 5 by p by n array giving the coverage at each simulation.
#' \code{COVERAGE[i, j, k]} is either 0 or 1 indicating weather the interval formed using
#' method i for rank j in simulation k covered its respective parameter}
#' \item{\code{WIDTH}}{ A 5 by p by n array giving the width of each interval in each simulation.}
#'  \item{\code{simnames}}{Name of each type of the 5 types of confidence intervals}
#' }
#'@export
cluster_sim <- function(beta, Sigma, err.sd=1, n.samp=100, n.rep=100, seed=NULL, parallel=FALSE){
  if(!is.null(seed)) set.seed(seed)
  
  lr_func <- function(data){
    y <- data[,1]
    X <- data[, -1]
    ests <- rcc:::many_lr(y, X, parallel=FALSE)
    df <- data.frame("estimate"=ests$beta_hat, "se"=ests$se_hat, "statistic"=ests$beta_hat/ests$se_hat)
    return(df)
  }
  
  
  rank_block <- function(stats, use.abs, blocks){
    p <- length(stats)
    b <- unique(blocks)
    N <- length(b)
    if(use.abs) stats <- abs(stats)
    rank <- rep(NA, p)
    top_ix_block <- t(sapply(b, FUN=function(blk){
      ix <- which(blocks==blk)
      ixmax <- which.max(stats[ix])
      return(c(ix[ixmax], max(stats[ix])))
    }))
    o <- order(top_ix_block[,2], decreasing=TRUE)
    j <- top_ix_block[order(top_ix_block[,2], decreasing=TRUE),1]
    rank <- match(1:p, j)
    return(list("order"=j, "rank"=rank))
  }
  
  K <- length(beta)
  stopifnot(length(Sigma) ==K)
  nK <- sapply(beta, FUN=length)
  stopifnot(all(nK == sapply(Sigma, FUN=nrow)))
  stopifnot(all(nK == sapply(Sigma, FUN=ncol)))
  p <- length(unlist(beta))
  blocks <- rep(1:K, nK)
  
  effects <- c()
  for(k in 1:K){
    effects <- c(effects, Sigma[[k]]%*%beta[[k]])
  }
  nms <- c("basic", "cw")
  ixs <- list(1:p, 1:K)
  simnames <- c("nonpar", "par", "wfb", "wfb2", "naive", "selInf1", "ash")
  simnames <- unlist(lapply(nms, FUN=function(nm){paste0(simnames, "_", nm)}))
 
  COVERAGE <- array(dim=c(length(simnames), p, n.rep))
  WIDTH <- array(dim=c(length(simnames),  p, n.rep))
  i <- 1
  while(i <= n.rep){
    cat(i, " ")
    #Generate data
    xs <- lapply(1:K, FUN=function(k){
      nk <- length(beta[[k]])
      mvrnorm(n=n.samp, mu = rep(0, nk), Sigma = Sigma[[k]])
    })
    X <- do.call(cbind, xs)
    y <- X%*% unlist(beta) + rnorm(n=n.samp, sd=err.sd)
    #Run linear regressions for each column of X
    data <- cbind(y, X)
    res.orig <- lr_func(data)
    #Ranking functions
    rk.basic <- rcc:::basic_rank(res.orig$statistic, use.abs=TRUE)
    rk.cw <- rank_block(stats=res.orig$statistic, use.abs=TRUE, blocks=blocks)
    orders <- list(rk.basic$order, rk.cw$order)
    cat("Got effect sizes.\n")
    
    #Non parametric bootstrap - basic
    ci.nonpar <- nonpar_bs_ci(data, analysis.func = lr_func,n.rep=1000, res.orig=res.orig, level = 0.9, parallel=parallel)$ci
    COVERAGE[which(simnames=="nonpar_basic"), ,i] <- (ci.nonpar[,1] < effects & effects < ci.nonpar[,2])[rk.basic$order]
    WIDTH[which(simnames=="nonpar_basic"), , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[rk.basic$order]
    
    
    #Non parametric bootstrap - cw
    ci.nonpar <- nonpar_bs_ci(data, analysis.func = lr_func,n.rep=1000, rank.func = rank_block,
                              res.orig=res.orig, level = 0.9, parallel=parallel, blocks=blocks)$ci
    COVERAGE[which(simnames=="nonpar_cw"), 1:K,i] <- (ci.nonpar[,1] < effects & effects < ci.nonpar[,2])[rk.cw$order]
    WIDTH[which(simnames=="nonpar_cw"), 1:K , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[rk.cw$order]
    
    
    #Parametric bootstrap - basic
    ci.par <- par_bs_ci(res.orig$estimate, res.orig$se,n.rep=1000, level=0.9)$ci
    COVERAGE[which(simnames=="par_basic"), ,i] <- (ci.par[,1] < effects & effects < ci.par[,2])[rk.basic$order]
    WIDTH[which(simnames=="par_basic"), , i] <- (ci.par[,2] -ci.par[,1])[rk.basic$order]
    
    #Parametric bootstrap - cw
    ci.par <- par_bs_ci(res.orig$estimate, res.orig$se,n.rep=1000, level=0.9, rank.func = rank_block, blocks=blocks)$ci
    COVERAGE[which(simnames=="par_cw"), 1:K ,i] <- (ci.par[,1] < effects & effects < ci.par[,2])[rk.cw$order]
    WIDTH[which(simnames=="par_cw"), 1:K, i] <- (ci.par[,2] -ci.par[,1])[rk.cw$order]
    
    #WFB CIs conditional on taking the top 10%
    t <- res.orig$estimate/res.orig$se
    ct <- abs(t)[rk.basic$order][floor(0.1*p) + 1]
    wfb <- lapply(t, FUN=function(x){
      if(abs(x) < ct) return(c(NA, NA))
      ci <- try(Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
      if(class(ci) == "try-error") return(c(NA, NA))
      return(ci)
    })
    ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=p)
    ci.wfb[,1] <- ci.wfb[,1]*res.orig$se
    ci.wfb[,2] <- ci.wfb[,2]*res.orig$se
    for(jj in 1:2){
      nm <- paste0("wfb_", nms[jj])
      COVERAGE[which(simnames==nm), ixs[[jj]] ,i] <- (ci.wfb[,1] < effects & effects < ci.wfb[,2])[orders[[jj]]]
      WIDTH[which(simnames==nm), ixs[[jj]], i] <- (ci.wfb[,2] -ci.wfb[,1])[orders[[jj]]]
    }
    
    #Naive ci
    ci.naive <- cbind(res.orig$estimate - res.orig$se*qt(0.95, df=p-1),
                      res.orig$estimate + res.orig$se*qt(0.95, df=p-1))
    for(jj in 1:2){
      nm <- paste0("naive_", nms[jj])
      COVERAGE[which(simnames==nm), ixs[[jj]] ,i] <- (ci.naive[,1] < effects & effects < ci.naive[,2])[orders[[jj]]]
      WIDTH[which(simnames==nm), ixs[[jj]], i] <- (ci.naive[,2] -ci.naive[,1])[orders[[jj]]]
    }    
    
    #Reid, Taylor, Tibshirani method (selectiveInference)
    M <- manyMeans(y=t, k=0.1*p, alpha=0.1, sigma=1)
    ci.rtt1 <- matrix(nrow=p, ncol=2)
    ci.rtt1[M$selected.set, ] <- M$ci
    ci.rtt1[,1] <- ci.rtt1[,1]*res.orig$se
    ci.rtt1[,2] <- ci.rtt1[,2]*res.orig$se
    for(jj in 1:2){
      nm <- paste0("selInf1_", nms[jj])
      COVERAGE[which(simnames==nm), ixs[[jj]] ,i] <- (ci.rtt1[,1] < effects & effects < ci.rtt1[,2])[orders[[jj]]]
      WIDTH[which(simnames==nm), ixs[[jj]], i] <- (ci.rtt1[,2] -ci.rtt1[,1])[orders[[jj]]]
    }
    
    #ashr credible intervals
    ash.res <- ash(betahat = res.orig$estimate, sebetahat = res.orig$se, mixcompdist = "normal")
    ci.ash <- ashci(ash.res, level=0.9, betaindex = 1:p, trace=FALSE)
    for(jj in 1:2){
      nm <- paste0("ash_", nms[jj])
      COVERAGE[which(simnames==nm), ixs[[jj]] ,i] <- (ci.ash[,1] < effects & effects < ci.ash[,2])[orders[[jj]]]
      WIDTH[which(simnames==nm), ixs[[jj]], i] <- (ci.ash[,2] -ci.ash[,1])[orders[[jj]]]
    }
    
    i <- i+1
    
  }
  
  cat("\n")
  return(list("COVERAGE"=COVERAGE,  "WIDTH"=WIDTH, "simnames"=simnames))
}

