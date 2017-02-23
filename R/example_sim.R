
#' Calculate expected coverage of several alternative intervals given a set of parameter values.
#' @description This function can be used to generate the example shown in figure 1.
#' @param theta Vector of parameter values (length p)
#' @param n Number of simulations to run.
#' @param use.abs Rank estimates by absolute value for bootstrapping methods.
#' @return A list with elements:
#' \describe{
#' \item{\code{COVERAGE}}{A 5 by p by n array giving the coverage at each simulation.
#' \code{COVERAGE[i, j, k]} is either 0 or 1 indicating weather the interval formed using
#' method i for rank j in simulation k covered its respective parameter}
#' \item{\code{WIDTH}}{ A 5 by p by n array giving the width of each interval in each simulation.}
#'  \item{\code{Z}}{A p by n matrix of all of the simulated data}
#'  \item{\code{simnames}}{Name of each type of the 5 types of confidence intervals}
#' }
#'@export
example_sim <- function(theta, n=200, use.abs=TRUE){
  p <- length(theta)
  Z <- matrix(nrow=p, ncol=n)
	simnames <- c("par", "wfb", "oracle", "naive", "selInf1", "ash")
	COVERAGE <- array(dim=c(length(simnames), p, n))
	WIDTH <- array(dim=c(length(simnames),  p, n))
	THETA <- array(dim=c(p, n))
  for(i in 1:n){
    cat(i, " ")
    z <- Z[,i] <- rnorm(n=p, mean=theta)
    if(use.abs) j <- order(abs(z), decreasing = TRUE)
      else j <- order(z, decreasing=TRUE)
    jinv <- match(1:p, j)
    THETA[,i] <- theta[j]
    #Parametric Bootstrap
    ci <- par_bs_ci(beta=z, theta=z, use.abs=use.abs, n=1000)[, c("ci.lower", "ci.upper")]
    COVERAGE[simnames == "par", ,i]<- ci[j,1] <= theta[j] & ci[j,2] >= theta[j]
		WIDTH[simnames=="par", , i] <- (ci[, 2] - ci[,1])[j]

    #Oracle
		ci <- par_bs_ci(beta=z, theta=theta, use.abs=use.abs, n=1000)[, c("ci.lower", "ci.upper")]
    COVERAGE[simnames == "oracle", ,i]<- ci[j,1] <= theta[j] & ci[j,2] >= theta[j]
		WIDTH[simnames=="oracle", , i] <- (ci[, 2] - ci[,1])[j]
    #Naive
    ci <- cbind(z-qnorm(0.95), z+qnorm(0.95))
    COVERAGE[simnames == "naive", ,i]<- ci[j,1] <= theta[j] & ci[j,2] >= theta[j]
		WIDTH[simnames=="naive", , i] <- (ci[, 2] - ci[,1])[j]

		#WFB CIs conditional on taking the top 10%
		#ct is the threshold on the absolute value of z -- no one tailed selection for wfb method
		if(use.abs) ct <- quantile(abs(z), probs=0.9)
		  else ct <- quantile(z, probs=0.9)
		wfb <- lapply(z, FUN=function(x){
			  if(abs(x) < ct) return(c(NA, NA))
			  ci <- try(rccSims::Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
			  if(class(ci) == "try-error") return(c(NA, NA)) #Sometimes WFB code produces errors
			  return(ci)
			 })
		ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=p)
    COVERAGE[simnames == "wfb", ,i]<- (ci.wfb[,1] <= theta & theta <= ci.wfb[,2])[j]
		WIDTH[simnames=="wfb", , i] <- (ci.wfb[, 2] - ci.wfb[,1])[j]

		#Reid, Taylor, Tibshirani method (selectiveInference) (known sigma)
		M <- manyMeans(y=z, k=0.1*p, alpha=0.1, sigma=1)
		ci.rtt1 <- matrix(nrow=p, ncol=2)
		ci.rtt1[M$selected.set, ] <- M$ci
		COVERAGE[simnames == "selInf1", ,i]<- (ci.rtt1[,1] <= theta & theta <= ci.rtt1[,2])[j]
		WIDTH[simnames=="selInf1", , i] <- (ci.rtt1[, 2] - ci.rtt1[,1])[j]
		#}

    #ASH
		ash.res <- ash(betahat = z, sebetahat = rep(1, p), mixcompdist = "normal")
		ci <- ashci(ash.res, level=0.9, betaindex = 1:p, trace=FALSE)
		COVERAGE[simnames == "ash", ,i]<- ci[j,1] <= theta[j] & ci[j,2] >= theta[j]
		WIDTH[simnames=="ash", , i] <- (ci[, 2] - ci[,1])[j]
  }
  cat("\n")
  return(list("Z"=Z, "COVERAGE"=COVERAGE, "WIDTH"=WIDTH,
              "simnames"=simnames, "THETA"=THETA))
}
