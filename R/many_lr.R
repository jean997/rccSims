#'@import tidyr
#'@import ggplot2
#'@import selectiveInference
#'@import ashr

#Run many linear regressions
many_lr <- function(y, X, parallel=FALSE){
  p <- dim(X)[2]
  n <- dim(X)[1]

  my_lr_func <- function(i){
    if(all(X[,i] == X[1, i])) return(c(0, 0))
    my.X <- cbind(rep(1, n), X[,i])
    Xsum <- sum(X[,i])
    Xsqsum <- sum(X[,i]^2)
    XtXinv <- solve(matrix( c(n, Xsum, Xsum, Xsqsum), nrow=2, byrow=T) )
    Xty <- matrix(c(sum(y), sum(y*X[,i])), nrow=2)
    beta_hat <- XtXinv%*%Xty
    res <- y - my.X %*% beta_hat
    s_hat <-  1/(n-2)*sum(res^2)
    V <- XtXinv%*%t(my.X) %*% diag(rep(s_hat, n)) %*% my.X %*% XtXinv
    return(c(beta_hat[2], sqrt(V[2,2])))
  }
  if(!parallel){
    Z <- lapply(1:p, FUN=function(i){
      my_lr_func(i)
    })
  }else{
    cores <- detectCores()-1
    cl <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(cl))
    Z <- parLapply(cl, 1:p, fun=function(i){
      my_lr_func(i)
    })
  }
  Z <- matrix(unlist(Z), nrow=2, ncol=p, byrow=FALSE)
  df <- data.frame("beta_hat"=Z[1,], "se_hat"=Z[2,])
  return(df)
}
