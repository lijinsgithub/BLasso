#' Bayesian Lasso by Variational Bayes + Dirichlet-Laplace Priors
#'
#' @description Extended Bayesian Lasso (Park \& Casella (2008)) with Dirichlet-Laplace Priors (Bhattacharya et al. 2015)
#'
#' @references Bhattacharya, Anirban, et al. "Dirichlet–Laplace priors for optimal shrinkage." Journal of the American Statistical Association 110.512 (2015): 1479-1490.
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
#' @examples ex1<-DLasso(x=data1[,-1], y=data1[,1]); ex1$beta; ex1$lambda; sum(ex1$p.value<0.05); #n of selected variables#
#'
#' @export
#'
DLasso<-function(x, y, print.it=FALSE){

  ##This is a function for Bayesian Lasso written for my dissertation ##
  ##The algoritm is implemented by Variational Bayes##

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

   psi <- 1
   alpha <- 1/2

  inv.D<-diag(as.vector(inv.tau2))
  A <- XtX +inv.D
  inv.A <- ginv(A)

  tol <-10^-2
  maxiter <- 10^2
  i<-0
  conv <-0


  L.comp <- function() {
    -((n+p-1)/2+1)*log(sigma2)-1/(2*sigma2)*sum(resid^2)-1/2*sum(log(tau2))-1/(2*sigma2)*sum(beta^2/(tau2))+sum(log(psi/2))-sum(psi/2*tau2)+(alpha-1)*I(alpha-1>0)*sum(log(psi)) - 1/2*sum(psi)
  }

  Approx <- function() {
    -((n+p-1)/2+1)*log(sigma21)-1/(2*sigma21)*sum(resid1^2)-1/2*sum(log(tau21))-1/(2*sigma21)*sum(beta1^2/(tau21))+sum(log(psi1/2))-sum(psi1/2*tau21)+(alpha1-1)*I(alpha1-1>0)*sum(log(psi1)) - 1/2*sum(psi1)
  }


  ELBO <- L.comp()
  diff1 <- ELBO

  if(print.it == TRUE){cat(" Iter. FE", alpha, ELBO , fill=T)}
  repeat {

    inv.D<-diag(as.vector(inv.tau2))
    A <- XtX +inv.D
    inv.A1 <- ginv(A)
    beta1 <- inv.A1%*%t(x)%*%y
    beta2 <- beta1^2
    xb <- x %*% beta1
    resid1 <- (y-xb)
    a <- (n+p+1)/2
    b <-t(resid1) %*% resid1/2 + 1/2* t(beta2) %*%inv.tau2

    sigma21 <-b/(a+1)

    inv.tau21 <- sqrt(psi^2*as.numeric(sigma21)/beta1^2)

    tau21 <- 1/inv.tau21 + 1/(psi^2)

    psi1 <- (alpha+1)/tau21

    psi1<-psi1/sum(psi1)

   s <- log(mean(psi1^2)) - mean(log(psi1^2))

    a.grad<- function(aa){
      return(log(aa) - digamma(aa))
    }
    a.hess <- function(aa){
      return(1/aa - trigamma(aa))
    }

     a0 <-(alpha+1)/2
     a.diff<-1

     while(a.diff > tol){
       a1 <- a0 - runif(1, 0, tol)*(a.grad(a0)-s)/a.hess(a0)
       a1 <- ifelse(a1<=0.5, a0+runif(1,0.5,1), a1)
       a.diff <- abs(a1-a0)
       a0 <-a1
   }
     alpha1 <-2*a1-1


    ELBO1 <- Approx()

    diff1 <- ELBO1 - ELBO

    if(print.it == TRUE){cat(" ", i," ", alpha1," "," ",ELBO1, fill=T)}

    if(i==maxiter) {
      conv<-0
      break}

    if(i>5 && diff1<tol) {
      conv <- 1
      break
    }

    i <- i + 1
    ELBO <- ELBO1
    beta<-beta1
    sigma2<-sigma21
    tau2<- tau21
    alpha <- alpha1
    inv.tau2 <- inv.tau21
    psi <- psi1
    resid <- resid1
    inv.A <- inv.A1
  }


  beta.sig<-sqrt(diag(inv.A)*sigma2)

  beta.t <-beta/beta.sig
  H <- x%*%inv.A%*%t(x)
  df.t <-  sum(diag(H))
  t.pval= round(apply(cbind(pt(beta.t, df.t), 1-pt(beta.t, df.t)),1, min), 3)

  list(beta=round(beta,3), beta.sig=beta.sig, tau2=tau2,  p.value=t.pval, sigma2=sigma2,alpha=alpha,  convergence = conv)
}
