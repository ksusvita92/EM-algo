#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to get EM-algo fro shape and scale
## Input: data, number of component
## Output: MoM estimator of shape and scale using EM-algo
## attr: class "gammaMix"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
EMgamma <- function(x, n.comp, alpha0=NULL, beta0=NULL, maxiter=1000, tol=1e-5, verbose=F) {
  if(n.comp>4) stop("We assume no more than 4 components in within a cluster.")
  if(is.null(alpha0)) {alpha <- (mean(x))/var(x)} #--- //use method of moment estimator
  else {alpha <- alpha0}
  if(is.null(beta0)) {beta <- mean(x)/var(x)} #--- //use method of moment estimator
  else {beta <- beta0}
  
  diff <- 1
  iter <- 1
  
  #--- //define component pdf:
  #---   Notice that we have to consider the discretization. So the pdf is f(x)g(x)
  #---   See reference for getting factor 2(1-x)
  f1_ <- function(x, a, b) (2-2*x)*dgamma(x, a, 1/b)
  f1upper_ <- function(x, a, b, r) (r+1-x)*dgamma(x, a, 1/b)
  f1lower_ <- function(x, a, b, r) (x-r+1)*dgamma(x, a, 1/b)
  
  f2_ <- function(x, a, b) (2-2*x)*dgamma(x, 2*a, 1/b)
  f2upper_ <- function(x, a, b, r) (r+1-x)*dgamma(x, 2*a, 1/b)
  f2lower_ <- function(x, a, b, r) (x-r+1)*dgamma(x, 2*a, 1/b)
  
  f3_ <- function(x, a, b) (2-2*x)*dgamma(x, 3*a, 1/b)
  f3upper_ <- function(x, a, b, r) (r+1-x)*dgamma(x, 3*a, 1/b)
  f3lower_ <- function(x, a, b, r) (x-r+1)*dgamma(x, 3*a, 1/b)
  
  f4_ <- function(x, a, b) (2-2*x)*dgamma(x, 4*a, 1/b)
  f4upper_ <- function(x, a, b, r) (r+1-x)*dgamma(x, 4*a, 1/b)
  f4lower_ <- function(x, a, b, r) (x-r+1)*dgamma(x, 4*a, 1/b)
  
  #--- //the discretization: P(T=t)=int_{t-1}^{t+1} f(x)g(x) dx
  p1 <- function(d, a, b) integrate(f=f1_, lower=d, upper=(d+1), a=a, b=b)
  p1upper <- function(d, a, b) integrate(f=f1upper_, lower=d, upper=(d+1), a=a, b=b, r=d)
  p1lower <- function(d, a, b) integrate(f=f1lower_, lower=max(d-1,0), upper=d, a=a, b=b, r=d)
  
  p2 <- function(d, a, b) integrate(f=f2_, lower=d, upper=(d+1), a=a, b=b)
  p2upper <- function(d, a, b) integrate(f=f2upper_, lower=d, upper=(d+1), a=a, b=b, r=d)
  p2lower <- function(d, a, b) integrate(f=f2lower_, lower=max(d-1,0), upper=d, a=a, b=b, r=d)
  
  p3 <- function(d, a, b) integrate(f=f3_, lower=d, upper=(d+1), a=a, b=b)
  p3upper <- function(d, a, b) integrate(f=f3upper_, lower=d, upper=(d+1), a=a, b=b, r=d)
  p3lower <- function(d, a, b) integrate(f=f3lower_, lower=max(d-1,0), upper=d, a=a, b=b, r=d)
  
  p4 <- function(d, a, b) integrate(f=f4_, lower=d, upper=(d+1), a=a, b=b)
  p4upper <- function(d, a, b) integrate(f=f4upper_, lower=d, upper=(d+1), a=a, b=b, r=d)
  p4lower <- function(d, a, b) integrate(f=f4lower_, lower=max(d-1,0), upper=d, a=a, b=b, r=d)
  
  tau1 <- numeric(length(x)); tau2 <- numeric(length(x))
  tau3 <- numeric(length(x)); tau4 <- numeric(length(x))
  while(diff>tol & iter<maxiter & !(is.na(diff))) {
    #--- //E-step: compute probability of interval belonging to a component
    for(id in 1:length(x)) {
      if(x[id]==0) {
        d1 <- p1(0, alpha, beta)[[1]] #--- //get the P(T=0)
        d2 <- p2(0, alpha, beta)[[1]]
        d3 <- p3(0, alpha, beta)[[1]]
        d4 <- p4(0, alpha, beta)[[1]]
      } else {
        d1 <- p1lower(x[id], alpha, beta)[[1]] + p1upper(x[id], alpha, beta)[[1]]
        d2 <- p2lower(x[id], alpha, beta)[[1]] + p2upper(x[id], alpha, beta)[[1]]
        d3 <- p3lower(x[id], alpha, beta)[[1]] + p3upper(x[id], alpha, beta)[[1]]
        d4 <- p4lower(x[id], alpha, beta)[[1]] + p4upper(x[id], alpha, beta)[[1]]
      }
      #--- //E-step: compute the relative probability of a data point to belong to 
      #---   one of the components
      if(n.comp==1) {all <- d1; tau1[id] <- d1/all}
      else if(n.comp==2) {all <- d1+d2; tau1[id] <- d1/all; tau2[id] <- d2/all}
      else if(n.comp==3) {all <- d1+d2+d3; tau1[id] <- d1/all; tau2[id] <- d2/all; tau3[id] <- d3/all}
      else if(n.comp==4) {all <- d1+d2+d3+d4; tau1[id] <- d1/all; tau2[id] <- d2/all; tau3[id] <- d3/all; tau4[id] <- d4/all}
    }
    
    #--- //E-step: get the weight omega
    if(n.comp==1) {w1 <- 1; w2 <- 0; w3 <- 0; w4 <- 0}
    else if(n.comp==2) {w1 <- mean(tau1); w2 <- mean(tau2); w3 <- 0; w4 <- 0}
    else if(n.comp==3) {w1 <- mean(tau1); w2 <- mean(tau2); w3 <- mean(tau3); w4 <- 0}
    else if(n.comp==4) {w1 <- mean(tau1); w2 <- mean(tau2); w3 <- mean(tau3); w4 <- mean(tau4)}

    
    #--- //M-step: update alpha and beta
    #---   we use weighted tau1 to compute: xbar = sum xi*pi
    update.xbar <- weighted.mean(x, tau1)
    cons_ <- sum(tau1)/((sum(tau1))^2-sum(tau1^2))
    update.var <- cons_*(sum(tau1*(x-update.xbar)^2))
    updatealpha <- update.xbar^2/update.var #--- //use MoM estimator
    updatebeta <- update.var/update.xbar #--- //use MoM estimator
    weights <- c(w1,w2,w3,w4)
    
    #--- //compute diff to check convergence
    diff <- sqrt((alpha-updatealpha)^2 + (beta-updatebeta)^2)
    alpha <- updatealpha
    beta <- updatebeta
    iter <- iter+1
    
    #--- //print updates
    if(verbose) {cat("Iteration", iter, ": alpha=", alpha, ", beta=", beta, ", diff=", diff, "\n")}
  }
  if(iter==maxiter) print("NOT CONVERGENT!")
  
  #--- //returns outputs
  names(weights) <- paste("w", 1:4, sep="")
  para <- c(alpha/n.comp, beta/w1)
  names(para) <- c("alpha", "beta")
  
  output <- list(weights=weights, parameter=para, diff=diff, niter=iter)
  class(output) <- "gammaMix"
  return(output)
}
EMgammaND <- function(x, n.comp, alpha0=NULL, beta0=NULL, maxiter=1000, tol=1e-5, verbose=F) {
  if(n.comp>4) stop("We assume no more than 4 components in within a cluster.")
  if(is.null(alpha0)) {alpha <- (mean(x))/var(x)} #--- //use method of moment estimator
  else {alpha <- alpha0}
  if(is.null(beta0)) {beta <- mean(x)/var(x)} #--- //use method of moment estimator
  else {beta <- beta0}
  
  diff <- 1
  iter <- 1
  
  tau1 <- numeric(length(x)); tau2 <- numeric(length(x))
  tau3 <- numeric(length(x)); tau4 <- numeric(length(x))
  
  #--- //if we have x==0
  f1_ <- function(x, a, b) (2-2*x)*dgamma(x, a, 1/b)
  f2_ <- function(x, a, b) (2-2*x)*dgamma(x, 2*a, 1/b)
  f3_ <- function(x, a, b) (2-2*x)*dgamma(x, 3*a, 1/b)
  f4_ <- function(x, a, b) (2-2*x)*dgamma(x, 4*a, 1/b)
  
  p1 <- function(d, a, b) integrate(f=f1_, lower=d, upper=(d+1), a=a, b=b)
  p2 <- function(d, a, b) integrate(f=f2_, lower=d, upper=(d+1), a=a, b=b)
  p3 <- function(d, a, b) integrate(f=f3_, lower=d, upper=(d+1), a=a, b=b)
  p4 <- function(d, a, b) integrate(f=f4_, lower=d, upper=(d+1), a=a, b=b)
  
  while(diff>tol & iter<maxiter & !(is.na(diff))) {
    #--- //E-step: compute probability of interval belonging to a component
    if(length(which(x==0))==0) {
      d1 <- pgamma(x, alpha, 1/beta)
      d2 <- pgamma(x, 2*alpha, 1/beta)
      d3 <- pgamma(x, 3*alpha, 1/beta)
      d4 <- pgamma(x, 4*alpha, 1/beta)
    } else {
      zerox_ <- which(x==0)
      nozerox <- which(x!=0)
      d1 <- numeric(length(x)); d2 <- numeric(length(x))
      d3 <- numeric(length(x)); d4 <- numeric(length(x))
      for(id in zerox_) {
        d1[id] <- p1(0, alpha, beta)[[1]]
        d2[id] <- p2(0, alpha, beta)[[1]]
        d3[id] <- p3(0, alpha, beta)[[1]]
        d4[id] <- p4(0, alpha, beta)[[1]]
      }
      d1[nozerox] <- pgamma(x[nozerox], alpha, 1/beta)
      d2[nozerox] <- pgamma(x[nozerox], 2*alpha, 1/beta)
      d3[nozerox] <- pgamma(x[nozerox], 3*alpha, 1/beta)
      d4[nozerox] <- pgamma(x[nozerox], 4*alpha, 1/beta)
    }
    
    if(n.comp==1) {all <- d1; tau1 <- d1/all}
    else if(n.comp==2) {all <- d1+d2; tau1 <- d1/all; tau2 <- d2/all}
    else if(n.comp==3) {all <- d1+d2+d3; tau1 <- d1/all; tau2 <- d2/all; tau3 <- d3/all}
    else if(n.comp==4) {all <- d1+d2+d3+d4; tau1 <- d1/all; tau2 <- d2/all; tau3 <- d3/all; tau4 <- d4/all}
    
    #--- //E-step: get the weight omega
    if(n.comp==1) {w1 <- 1; w2 <- 0; w3 <- 0; w4 <- 0}
    else if(n.comp==2) {w1 <- mean(tau1); w2 <- mean(tau2); w3 <- 0; w4 <- 0}
    else if(n.comp==3) {w1 <- mean(tau1); w2 <- mean(tau2); w3 <- mean(tau3); w4 <- 0}
    else if(n.comp==4) {w1 <- mean(tau1); w2 <- mean(tau2); w3 <- mean(tau3); w4 <- mean(tau4)}
    
    #--- //M-step: update alpha and beta
    #---   we use weighted tau1 to compute: xbar = sum xi*pi
    update.xbar <- weighted.mean(x, tau1)
    cons_ <- sum(tau1)/((sum(tau1))^2-sum(tau1^2))
    update.var <- cons_*(sum(tau1*(x-update.xbar)^2))
    updatealpha <- (update.xbar^2/update.var) #--- //use MoM estimator
    updatebeta <- (update.var/update.xbar) #--- //use MoM estimator
    weights <- c(w1,w2,w3,w4)
    
    #--- //compute diff to check convergence
    diff <- sqrt((alpha-updatealpha)^2 + (beta-updatebeta)^2)
    alpha <- updatealpha
    beta <- updatebeta
    iter <- iter+1
    
    #--- //print updates
    if(verbose) {cat("Iteration", iter, ": alpha=", alpha, ", beta=", beta, ", diff=", diff, "\n")}
  }
  if(iter==maxiter) print("NOT CONVERGENT!")
  
  #--- //returns outputs
  names(weights) <- paste("w", 1:4, sep="")
  para <- c(alpha/n.comp, beta/w1) #--- //this is scaling: since we use component 1 for estimation, 
  # notice that w1f1~Gamma(alpha, beta/w1) and for alpha/4 is bc there are 4 components
  names(para) <- c("alpha", "beta")
  
  output <- list(weights=weights, parameter=para, diff=diff, niter=iter)
  class(output) <- "gammaMix"
  return(output)
}
