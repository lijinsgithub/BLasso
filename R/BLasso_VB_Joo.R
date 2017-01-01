#' Bayesian Lasso by Variational Bayes
#'
#' @description Fit Bayesian Lasso (Park \& Casella (2008)) by variational Bayes
#'
#' @references Park, Trevor, and George Casella. "The bayesian lasso." Journal of the American Statistical Association 103.482 (2008): 681-686.
#' @references Joo, Lijin. "Bayesian Lasso: An Extension for Genome-wide Assoication Study." New York University, 2017
#'
#' @param x: predictor variables (numertic only)
#' @param y: outcome (numertic only)
#' @param print.it = TRUE/FALSE (default: FALSE, suppressing to print the number of iterations)
#'
#' @return beta
#' @return beta.sig: standard deviation of beta
#' @retrun tau2: local scale
#' @return p.value: p-value for t-test for beta (for variable selection)
#' @return sigma2
#' @return lambda: penalty, or global scale
#' @return convergence: 1/0 converged if convergence = 1
#'
#' @examples ex1<-VBLasso(x=data1[,-1], y=data1[,1]); ex1$beta; ex1$lambda; sum(ex1$p.value<0.05); #n of selected variables#
#'
#' @export
#'
VBLasso<-function(x, y, print.it = FALSE){

  ##as of 12.31##
  ##for fitting narrow regions, some adaptations are necessary##
  require(mvtnorm)
  require(MASS)
  require(pscl)
  require(statmod)

  n <- nrow(x)
  p <- ncol(x)
  x <- as.matrix(x, ncol=p, byrow=TRUE)

  meanx <- apply(x, 2, mean)
  x <- scale(x, meanx, FALSE)
  mu <- mean(y)
  y <- drop(y - mu)
  XtX <- t(x) %*% x
  xy <- t(x) %*% y

  #initial values#
  beta <- drop(backsolve(XtX + diag(nrow=p), xy))
  resid <- drop(y - x %*% beta)
  sigma2 <- drop((t(resid) %*% resid) / n)
  tau2 <- 1 / (beta * beta)
  inv.tau2 <- 1/tau2
  lambda2 <- p * sigma2 / sum(beta^2)


  tol <-10^-2
  maxiter <- 1000
  i<-0
  conv <-0

  L.comp <- function(){
    -((n+p-1)/2+1)*log(sigma2)-1/(2*sigma2*(sum(resid^2)))-1/2*sum(log(tau2))-1/2*sum(beta^2/(tau2))/sigma2+p*log(lambda2)-lambda2/2*sum(tau2)
  }

  Approx <- function() {
    -((n+p-1)/2+1)*log(sigma21)-1/(2*sigma21*(sum(resid1^2)))-1/2*sum(log(tau21))-1/2*sum(beta1^2/(tau21))/sigma21+p*log(lambda21)-lambda21/2*sum(tau21)
    }

  ELBO <- L.comp()
  diff1 <- ELBO
  lam.list<- c(lambda2)
  if(print.it == TRUE){cat(" Iter", "lam", "ELBO" , fill=T)}
  repeat {

    inv.D<-diag(as.vector(inv.tau2))
    A <- XtX +inv.D
    inv.A <- ginv(A)
    beta1 <- inv.A%*%t(x)%*%y
    beta2 <- beta1^2
    xb <- x %*% beta1
    resid1 <- (y-xb)
    a <- (n+p+1)/2
    b <-t(resid1) %*% resid1/2 + 1/2* t(beta2) %*%inv.tau2

    sigma21 <-b/(a+1)

    inv.tau21 <- sqrt(lambda2*as.numeric(sigma21)/beta1^2)

    tau21 <- 1/inv.tau21 + 1/lambda2


    lam.gradient <-function(lam, tt) {
      p/lam-1/2*sum(tau21)
    }
    lam.hessian <- function(lam){
      -p/(lam)^2
    }

    lam0 <-ifelse(lambda2>10, runif(1, 0, 10),lambda2)


    lam.diff<-1

    while(lam.diff>tol){
      lam1 <- lam0 -runif(1, 0,tol)*lam.gradient(lam0, tau21)/lam.hessian(lam0)
      lam.diff <- abs(lam1 -lam0)
      lam0<- lam1
    }

    lambda21 <- lam1

    ELBO1 <- Approx()
    diff1 <- ELBO1 - ELBO

    if(print.it == TRUE){cat(" ", i," ", sqrt(lambda21)," "," ",ELBO1, fill=T)}

    if(i==maxiter) {
      conv<-0
      break}

    if(i>5 && diff1< tol) {
      conv <- 1
      break
    }

    i <- i + 1
      ELBO <- ELBO1
      beta<-beta1
      sigma2<-sigma21
      tau2<- tau21
      lambda2 <- lambda21
      inv.tau2 <- inv.tau21
      resid <- resid1
  }

  beta.sig<-sqrt(diag(inv.A)*sigma2)

  beta.t <-beta/beta.sig
  H <- x%*%inv.A%*%t(x)
  df.t <-  sum(diag(H))
  t.pval= round(apply(cbind(pt(beta.t, df.t), 1-pt(beta.t, df.t)),1, min), 3)


  list(beta=round(beta,3), beta.sig=beta.sig, tau2=tau2, p.value=t.pval, sigma2=sigma2,lambda=sqrt(lambda2), convergence = conv)
}
