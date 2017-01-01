#' Bayesian Lasso by Gibbs Sampler
#' @description Provide two options for the estimation of penalty parameter "lambda"
#'
#' @references Park, Trevor, and George Casella. "The bayesian lasso." Journal of the American Statistical Association 103.482 (2008): 681-686.
#'
#' @param x: predictor variables (numertic only)
#' @param y: outcome (numertic only)
#' @param n.max: n of interations (default: 10000)
#' @param EB: TRUE/FALSE (default: TRUE, estimating lambda by empircal bayes)
#' @param a, b: specify these for a hyper-Gamma prior for lambda^2 if EB = FALSE
#' @param print.it = TRUE/FALSE (default: FALSE, suppressing to print the number of iterations)
#'
#' @return beta
#' @return beta.95q: 95 \% (posterior) CI of beta
#' @return beta.sig: standard deviation of beta
#' @retrun tau2: local scale
#' @return tau2.95q: 95 \% CI of tau2
#' @retrun sigma2
#' @return sigma.95q: 95 \% CI of sigma2
#' @return lambda: penalty, or global scale
#' @return lambda.95q: 95% CI of lambda
#' @examples ex1<-BLasso(x=data1[,-1], y=data1[,1]); ex1$beta; ex1$lambda; sum(ex1$beta.95q[1,]*ex1$beta.95q[2,]>0); #n of selected variables#
#' @examples ex2 <-BLasso(x=data1[,-1], y=data1[,1], EB=FALSE, a=1.78, b=1);
#' @export
#'
BLasso = function(x, y, n.max = 10000, EB = TRUE, a=1, b=1, print.it=FALSE) {


  require("mvtnorm") #sampling from a multivarate normal#
  require("statmod") #sampling from inverse gaussian #
  require("pscl")    #samping from inverse gamma#
  require("Matrix")  #for fast matrix inversion#
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  #scaling data#
  n <- length(y)
  meanx <- drop( rep(1, n) %*% x)/n
  x <- scale(x, meanx, FALSE)
  mu <- mean(y)
  y <- drop(y - mu)

  XtX <- t(x) %*% x
  xy <- t(x) %*% y

  beta.sim <- matrix(0, n.max, p)
  sigma2.sim <- rep(0, n.max)
  tau2.sim <- matrix(0, n.max, p)
  lambda.sim <- rep(0, n.max)

  #initial values#
  beta <- drop(backsolve(XtX + diag(nrow=p), xy))
  resid <- drop(y - x %*% beta)
  sigma2 <- drop((t(resid) %*% resid) / n)
  invtau2 <- 1 / (beta * beta)
  lambda <- p * sqrt(sigma2) / sum(beta^2)

  if(EB==TRUE){

    a <-0;
    b <-0;
    tol <-10^-3
    Q<-function(){
      p*log(lambda^2)-lambda^2/2*sum(1/invtau2)
    }
      L0 <- Q()
  }

  iter <- 1
  while (iter < n.max) {
    if (iter %% round(n.max/4) == 0) {
   if(print.it==TRUE) cat('Iter:', iter, "", "\r")
    }

    # update beta
    invD <- diag(invtau2)
    invA <- ginv(XtX + invD)
    beta.m <- invA %*% xy
    Sigma <- as.numeric(sigma2) * invA
    if(det(Sigma)<0.001){
      beta <- drop(rmvnorm(1, beta.m, Sigma, "svd"))
    } else{
      beta <- drop(rmvnorm(1, beta.m, Sigma, "chol"))
    }

        beta.sim[iter,] <- beta

    # update sigma2
    sig.a <- (n+p-1)/2+a
    resid <- drop(y - x %*% beta)
    sig.b <- (t(resid) %*% resid + t(beta) %*% invD %*% beta)/2+b
    sigma2 <- rigamma(1, alpha=sig.a, beta=sig.b) #change inv-gamma function#
    sigma2.sim[iter] <- sigma2

    # update tau2
    mu.t <- sqrt(lambda^2 * sigma2 / beta^2)
    lambda.t <- lambda^2
    invtau2 <- rinvgauss(p, mean=mu.t, shape=lambda.t)
    tau2.sim[iter, ] <- 1/invtau2

    # update lambda

   if(EB==TRUE){

     lambda1 <- sqrt((2*p)/sum(1/invtau2))
     L1 <-Q()
     if(abs(L0-L1)>tol){
       lambda <- lambda1
       L0 <- L1
     } else {
       lambda <- lambda
       L0 <- L0
     }
     lambda.sim[iter]<-lambda

   } else{

     sh <- p+a
     sc <- sum(1/invtau2)/2+b
     lambda<-sqrt(rgamma(1, shape=sh, rate=sc))

     lambda.sim[iter]<-lambda
   }


    iter <- iter + 1

  }

  colnames(beta.sim) <- colnames(x)
  colnames(tau2.sim) <- colnames(x)

  list(beta=round(apply(beta.sim[seq(round(n.max/2), n.max),],2, median),3),
       beta.95q = round(apply(beta.sim[seq(round(n.max/2), n.max),], 2, function(x){quantile(x,c(0.025, 0.975))}),3),
       beta.sig = round(apply(beta.sim[seq(round(n.max/2), n.max),], 2, sd),3),
       tau2=round(apply(tau2.sim[seq(round(n.max/2), n.max),],2,  median),3),
       tau2.95q = round(apply(tau2.sim[seq(round(n.max/2), n.max),], 2, function(x){quantile(x,c(0.025, 0.975))}),3),
       sigma2=round(median(sigma2.sim[seq(round(n.max/2), n.max)]),3),
       sigma2.95q = round(quantile(sigma2.sim[seq(round(n.max/2), n.max)], c(0.025, 0.975)),3),
       lambda=round(median(lambda.sim[seq(round(n.max/2), n.max)]),3),
       lambda.95q = round(quantile(lambda.sim[seq(round(n.max/2), n.max)], c(0.025, 0.975)),3))
}





