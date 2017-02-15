# Functions for calculating conditional acceptance regions of level alpha,
# and (1-alpha)100% conditional confidence intervals, 
# for normally distributed estimator with variance 1
# conditioned upon exceeding in absolute value a threshold ct.
# Escorting the paper "Selection Adjusted Confidence Intervals 
#  with More Power to Determine the Sign"
# By Asaf Weinstein, William Fithian, and Yoav Benjamini JASA (2012)


# Conditional Shortest Acceptance Region

Shortest.AR <- function(theta,ct,alpha){
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)	
#compute useful quantities
f <- function(theta) ( pnorm(ct + theta) - pnorm(ct - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0,ct + qnorm(1 - alpha))) $ root
f <- function(theta) 2*pnorm(theta - ct) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,ct + qnorm(1 - alpha/2) ) ) $ root      # theta2 - ct cannot be greater than qnorm(1 - alpha/2)!

#compute ends of the AR
is.neg <- 0
  if (theta < 0) is.neg <- 1
theta <- abs(theta)

  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -ct
  lr <- ct
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }

  if (0 < theta && theta < theta1) {
  ll <- theta - qnorm( 1 - alpha/2 *  Q(theta) )
  ul <- -ct
  lr <- ct
  ur <- theta + qnorm( 1 - alpha/2 *  Q(theta) )
    }

  if (theta1 <= theta && theta < theta2) {
  ll <- -ct
  ul <- -ct
  lr <- ct
  ur <- theta + qnorm( pnorm( ct - theta ) + (1 - alpha) *  Q(theta) )
    }

  if (theta2 <= theta) {
  ll <- -ct
  ul <- -ct
  lr <- theta - qnorm( .5 * ( 1 + (1 - alpha) *  Q(theta) ) )
  ur <- theta + qnorm( .5 * ( 1 + (1 - alpha) *  Q(theta) ) )
    }
l <- (ul - ll) + (ur - lr)  #compute the length of the AR
if(ll==-ct) {
	ll <- NA
	ul <- NA
}
A <- c(ll,ul,lr,ur)
if (is.neg == 1) A <- c(-ur,-lr,-ul,-ll)
v <- list(A = A,l = l)
return(v)
}



# Conditional Shortest Confidence Interval

Shortest.CI <- function(x,ct,alpha) {
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)	
#compute useful quantities
f <- function(theta) ( pnorm(ct + theta) - pnorm(ct - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, ct + qnorm(1 - alpha))) $ root
f <- function(theta) 2*pnorm(theta - ct) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,ct + qnorm(1 - alpha/2) ) ) $ root      # theta2 - ct cannot be greater than qnorm(1 - alpha/2)!
x0 <- qnorm( 1 - .5*alpha * Q(0) )  # not used here, but leaving it in
x1 <-  theta1 + qnorm ( pnorm(ct - theta1) + (1 - alpha)*( 1 - pnorm(ct + theta1) + 1 - pnorm(ct - theta1) ) )
x2 <- theta2 + (theta2 - ct)
f <- function(theta) 1 + ( -dnorm(ct - theta) + (1 - alpha) * ( -dnorm(ct + theta) + dnorm(ct - theta) ) ) / ( dnorm( qnorm( pnorm(ct - theta) + (1 - alpha) * (1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)) ) ) )
R <- uniroot(f, c(theta1, theta2)) $ root

#obtain CI ends
is.neg <- 0
  if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI
if (ct < x && x < x1) {
  f <- function(theta) 2 * (1 - pnorm(x - theta)) - alpha *  Q(theta)
  lower <- uniroot(f, c(-theta1 - .1, theta1)) $ root
  }
if (x1 < x && x < x2) {
  f <- function(theta) pnorm( x - theta ) - pnorm( ct - theta ) - (1 - alpha) *  Q(theta)
  lower <- uniroot(f, c(R, theta2)) $ root
  }
if (x > x2) {
  f <- function(theta) 2 * pnorm( x - theta ) - 1 - (1 - alpha) *  Q(theta)
  lower <- uniroot(f, c(theta2, x)) $ root
  }
#obtain upper end of CI
f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) *  Q(theta)
upper <- uniroot( f, c( theta2, x + 2 * qnorm(1 - alpha/2) ) ) $ root # upper is < x + 1.96, but for root-finding issues we put 2 * 1.96
CI <- list(lower = lower, upper = upper)
if (is.neg == 1) CI <- list(lower = -upper, upper = -lower)
return(CI)
}

# Conditional Modified Pratt Acceptance Region

# r is the allowed expansion of the length relative the to the conventional CI

MP.AR <- function(theta,r=1.2,ct,alpha){
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)
#compute useful quantities
f <- function(theta) pnorm( ct + r * Shortest.AR(theta,ct,alpha)$l - theta ) - pnorm( ct - theta ) - (1 - alpha) *  Q(theta)
thetatilde1 <- uniroot(f,c(0,theta1)) $ root
thetatilde2 <- uniroot(f,c(theta1, r * 2 * qnorm(1 - alpha / 2)))$root

#compute ends of the AR
is.neg <- 0
  if (theta < 0) is.neg <- 1
theta <- abs(theta)

  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -ct
  lr <- ct
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
  if (0 < theta && theta < thetatilde1) {
l <- Shortest.AR(theta,ct,alpha)$l
f <- function(x) 1 - pnorm(theta - x) + 1 - pnorm(ct + r * l - (-ct - x) - theta) - alpha *  Q(theta)
atilde1 <- uniroot(f,c(theta - qnorm(1 - alpha/2 *  Q(theta)), -ct))$root
  ll <- atilde1
  ul <- -ct
  lr <- ct
  ur <- ct + r * l - (-ct - atilde1)
    }
  if (thetatilde1 <= theta && theta < thetatilde2) {
l <- Shortest.AR(theta,ct,alpha)$l
f <- function(x) pnorm(x + r * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
atilde2 <- uniroot(f,c(ct, theta - qnorm( (1 - alpha) *  Q(theta) )))$root
  ll <- NA
  ul <- NA
  lr <- atilde2
  ur <- atilde2 + r * l
    }
  if (thetatilde2 <= theta) {
l <- Shortest.AR(theta,ct,alpha)$l
f <- function(x) pnorm(x + r * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
atilde2 <- uniroot(f,c(theta - r * l /2, theta - qnorm((1 - alpha) *  Q(theta))))$root
  ll <- NA
  ul <- NA
  lr <- atilde2
  ur <- atilde2 + r * l
    }
A <- c(ll,ul,lr,ur)
if (is.neg == 1) A <- c(-ur,-lr,-ul,-ll)
return(A)
}



# Conditional Modified Pratt Confidence Interval

# r is the allowed expansion of the length of acceptance region 
# relative the to the conventional CI

MP.CI <- function(x,r=1.2,ct,alpha){
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)	
#compute useful quantities
f <- function(theta) pnorm( ct + r * Shortest.AR(theta,ct,alpha)$l - theta ) - pnorm( ct - theta ) - (1 - alpha) *  Q(theta)
thetatilde1 <- uniroot(f,c(0,theta1)) $ root
thetatilde2 <- uniroot(f,c(theta1, r * 2 * qnorm(1 - alpha / 2)))$root
zalphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
lzero <- Shortest.AR(0,ct,alpha)$l
f <- function(x) 1 - pnorm(x) + 1 - pnorm(2 * ct + r * lzero - x) - 2 * alpha * (1 - pnorm(ct))
xtilde1 <- uniroot(f,c(ct,zalphahalf))$root
xtilde2 <- uniroot(f,c(zalphahalf,ct + r * 2 * qnorm(1 - alpha / 2)))$root
ltilde1 <- Shortest.AR(thetatilde1,ct,alpha)$l

#obtain CI ends
is.neg <- 0
  if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI
if (ct < x && x < xtilde1) {
  f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * ct + r * Shortest.AR(theta,ct,alpha)$l) - alpha *  Q(theta)
  lower <- uniroot(f,c(-thetatilde1 - 1e-3,0))$ root
  }
if (xtilde1 <= x && x < zalphahalf) {
  lower <- 0
  }
if (zalphahalf <= x && x < xtilde2) {
  lower <- 0
  }
if (xtilde2 <= x && x < ct + r * ltilde1 ) {
  f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * ct + r * Shortest.AR(theta,ct,alpha)$l) - alpha *  Q(theta)
  lower <- uniroot(f,c(0,thetatilde1))$ root
  }
if (x > ct + r * ltilde1 ) {
  f <- function(theta) pnorm(x - theta) - pnorm(x - r * Shortest.AR(theta,ct,alpha)$l - theta) - (1 - alpha) *  Q(theta)
  m <- optimize(f, c(thetatilde1,x), maximum=T)$maximum
  lower <- uniroot(f, c(thetatilde1,m))$root
  }
#obtain upper end of CI
f <- function(theta) pnorm(x + r * Shortest.AR(theta,ct,alpha)$l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
m <- optimize(f, c(x, x + r * 2 * qnorm(1 - alpha / 2)), maximum = T)$maximum
upper <- uniroot(f,c(0,m))$root
CI <- list(lower = lower, upper = upper)
if (is.neg == 1) CI <- list(lower = -upper, upper = -lower)
return(CI)
}


# Conditional Quasi-Conventinal Acceptance regions

# lambda governs the ballance between more power to determine the sign (smaller lambda) 
# and shorter length of acceptance region (larger lambda)

QC.AR <- function(theta,lambda=0.4,ct,alpha){
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)
#compute useful quantities
f <- function(theta) ( pnorm(ct + theta) - pnorm(ct - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, ct + qnorm(1 - alpha))) $ root
f <- function(theta) 1 - pnorm(ct - theta) - (1 - alpha) * Q(theta)
thetastar <- uniroot(f,c(0,ct + qnorm(1 - alpha)))$root
f <- function(theta) 1 + lambda * ( 1 - dnorm(ct + theta) / dnorm(qnorm(2 - pnorm(ct + theta) - alpha * Q(theta))) )
b <- thetastar
while(is.na(f(b))==T) {b <- b + 1e-5}
thetaprime1 <- uniroot(f, c(b, theta1))$root

#compute ends of the AR
is.neg <- 0
  if (theta < 0) is.neg <- 1
theta <- abs(theta)

  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -ct
  lr <- ct
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
  if (0 < theta && theta < thetaprime1) {
f <- function(d) 1 + lambda * ( 1 - dnorm(ct + d + theta) / dnorm(qnorm(2 - pnorm(ct + d + theta) - alpha * Q(theta))) )
dunderbar <- max(-ct - theta + qnorm(1 - alpha * Q(theta)), 0)
dbar <- -ct - Shortest.AR(theta,ct,alpha)$A[1]
dstar <- uniroot(f,c(dunderbar+1e-5, dbar))$root
aprime1 <- theta - ct + qnorm( 1 - pnorm(ct + dstar + theta) + 1 - alpha * Q(theta) )
  ll <- -ct - dstar
  ul <- -ct
  lr <- ct
  ur <- ct + aprime1
    }
  if (thetaprime1 <= theta && theta < theta1) {
aprime2 <- theta - ct + qnorm( pnorm(ct - theta) + (1 - alpha) * Q(theta) )
  ll <- NA
  ul <- NA
  lr <- ct
  ur <- ct + aprime2
    }
  if (theta > theta1) {
A <- Shortest.AR(theta,ct,alpha)$A
  ll <- A[1]
  ul <- A[2]
  lr <- A[3]
  ur <- A[4]
    }
A <- c(ll, ul, lr, ur)
if (is.neg == 1) A <- c(-ur,-lr,-ul,-ll)
return(A)
}



#Conditional Quasi-Conventinal Confidence Interval 

# lambda governs the ballance between more power to determine the sign (smaller lambda) 
# and shorter length of the confidence interval (larger lambda)


QC.CI <- function(x,lambda=0.4,ct,alpha){
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)
f <- function(theta) ( pnorm(ct + theta) - pnorm(ct - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, ct + qnorm(1 - alpha))) $ root
f <- function(theta) 1 - pnorm(ct - theta) - (1 - alpha) * Q(theta)
thetastar <- uniroot(f,c(0,ct + qnorm(1 - alpha)))$root
f <- function(theta) 1 + lambda * ( 1 - dnorm(ct + theta) / dnorm(qnorm(2 - pnorm(ct + theta) - alpha * Q(theta))) )
b <- thetastar
while(is.na(f(b))==T) {b <- b + 1e-5}
thetaprime1 <- uniroot(f, c(b, theta1))$root

dmin <- function(theta) max( -ct - theta + qnorm(1 - alpha * Q(theta)), 0 )
dmax <- function(theta) -ct - Shortest.AR(theta,ct,alpha)$A[1]
#compute useful quantities
zalphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
dminzero <- dmin(0)
dmaxzero <- dmax(0)
f <- function(d) 1 + lambda * ( 1 - dnorm(ct + d) / dnorm(qnorm(2 - pnorm(ct + d) - alpha * Q(0))) )
xprime1 <- uniroot(f, c(dminzero, dmaxzero))$root + ct
xprime2 <- qnorm( 2 - pnorm(xprime1) - alpha * Q(0) )
xprime3 <- thetaprime1 + qnorm( pnorm(ct - thetaprime1) + (1 - alpha) * Q(thetaprime1) )

#obtain CI ends
is.neg <- 0
if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI        
#NOTE: numerical problem may arise in obtaining lower end when |x-ct| extremely small 
# (this problem was not encountered in our use of the fuction)
if (ct < x && x < xprime1) {
d <- x - ct
f <- function(theta) 1 + lambda * ( 1 - dnorm(ct + d + theta) / dnorm(qnorm(2 - pnorm(ct + d + theta) - alpha * Q(theta))) )
g <- function(theta) dmin(theta) - (x - ct)
  if (x - ct > dmin(0)) {thetamin <- 0} else {
    g <- function(theta) dmin(theta) - (x - ct)
    thetamin <- uniroot(g, c(0, thetastar + 1e-3))$root
}
leftend <- thetamin + 1e-5
while (is.na(suppressWarnings(f(leftend)))) leftend <- leftend + 1e-5 * 5
lower <- uniroot(f,c(leftend, thetaprime1 + .1)) $ root       
#NOTE: for very small x, this numerical computation might find the lower bound 
# to be bigger than thetaprime1, 
# even though the true root must be smaller than thetaprime1 for such x; 
# However, the discrepancy is very small and the constrain is applied.
lower <- - lower
  }
if (xprime1 <= x && x < zalphahalf) {
  lower <- 0    # 0 is included in CI
  }
if (zalphahalf <= x && x < xprime2) {
  lower <- 0    # 0 is NOT included in CI
  }
if (xprime2 <= x && x < xprime3) {
g <- function(theta) -ct - theta + qnorm( 2 - pnorm(x - theta) - alpha * Q(theta) )
f <- function(theta) 1 + lambda * ( 1 - dnorm(ct + g(theta) + theta) / dnorm(qnorm(2 - pnorm(ct + g(theta) + theta) - alpha * Q(theta))) )
h <- function(theta) 1 - pnorm(x - theta) - alpha * Q(theta)
thetamax <- uniroot(h,c(0,x))$root
lower <- uniroot(f,c(0 - 1e-1, thetamax - 1e-1)) $ root
  }
if (x >= xprime3 ) {
  lower <- Shortest.CI(x,ct,alpha)$lower
  }
#obtain upper end of CI
f <- function(theta) 2*pnorm(theta - ct) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,ct + qnorm(1 - alpha/2) ) ) $ root      # theta2 - ct cannot be greater than qnorm(1 - alpha/2)!
f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) * Q(theta)
upper <- uniroot( f, c( theta2, x + 1.1 * qnorm(1 - alpha/2) ) ) $ root # upper is < x + 1.96, but for root-finding issues I put 2 * 1.96
CI <- list(lower = lower, upper = upper)
if (is.neg == 1) CI <- list(lower = -upper, upper = -lower)
return(CI)
}



# Maximum length increase of the Quasi-Conventional CI 
# as a function of lambda, for a given threshold ct and alpha.

QC.maxlength <- function(lambda,ct,alpha){
Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)
dmin <- function(theta) max( -ct - theta + qnorm(1 - alpha * Q(theta)), 0 )
dmax <- function(theta) -ct - Shortest.AR(theta,ct,alpha)$A[1]
f <- function(d) 1 + lambda * ( 1 - dnorm(ct + d) / dnorm(qnorm(2 - pnorm(ct + d) - alpha * Q(0))) )
dminzero <- dmin(0)
dmaxzero <- dmax(0)
while (is.na(suppressWarnings(f(dminzero)))) dminzero <- dminzero + 1e-12  # A numerical problem in evaluating f() at dminzero (or a value extremely close to it)
xprime1 <- uniroot(f, c(dminzero, dmaxzero))$root + ct

xprime2 <- qnorm( 2 - pnorm(xprime1) - alpha * Q(0) )
a <- QC.CI(xprime2,lambda,ct,alpha) $ lower
b <- QC.CI(xprime2,lambda,ct,alpha) $ upper
l <- b - a
return(l)
}